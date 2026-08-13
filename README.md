```
  _____                _          __          ___                  _
 |  __ \              | |         \ \        / (_)                | |              /\
 | |__) |_ _ _ __ __ _| | ___   __ \ \  /\  / / _ ______ _ _ __ __| |             /  \
 |  ___/ _` | '__/ `_ | |/ _ \ / _` \ \/  \/ / | |_  / `_ | '__/ _` |            |    |
 | |  | (_| | | | (_| | | (_) | (_| |\  /\  /  | |/ / (_| | | | (_| |          --:'''':--
 |_|   \__,_|_|  \__,_|_|\___/ \__, | \/  \/   |_/___\__,_|_|  \__,_|            :'_' :
                               |___/                                             _:"":\___
                                                                 ' '      ____.' :::     '._
                                                                . *=====<<=)           \    :
                                                                 .  '      '-'-'\_      /'._.'
                                                                                  \====:_ ""
                 _      ___   _____    ___                                       .'     \\
                /_\    / __| |_   _|  / __|                                     :       :
               / _ \  | (__    | |   | (_ |                                    /   :    \
              /_/ \_\  \___|   |_|    \___|                                   :   .      '.
                                                              ,. _            :  : :      :
           .   .-. .-.   .-. .-.   .-. .-.   .             '-' _  ).          :__:-:__.;--'
           |\ /|||\|||\ /|||\|||\ /|||\|||\ /|           (   _|  _  )        '-'   '-'
           ||\|||/ \|||\|||/ \|||\|||/ \|||\||        ( -  _| |_|   -_
           -~ `-~   `-~ `-`   `-~ `-`   `-~ `-       (   _| |_  |_    )
                                                     '-   |_         -
```

**ParalogWizard** detects and separates paralogs in targeted HybSeq / exon-capture data, then builds customized references and orthologous alignments.

Each cast step **overwrites** its outputs. Pass `--debug` on any subcommand for DEBUG lines in the shared `*.log` file and on the console.



---



# Dependencies

**Python** (3.8+ recommended) with:

- Biopython, NumPy, SciPy, scikit-learn, Matplotlib, pandas

**External tools** (must be on `PATH`):


| Step                       | Tools                                          |
| -------------------------- | ---------------------------------------------- |
| `cast_assemble`            | BWA, samtools, SPAdes                          |
| `cast_retrieve`            | makeblastdb, blastn                            |
| `cast_analyze`             | MAFFT, FastTree (or FastTreeMP / VeryFastTree) |
| `cast_separate`            | BLAT, MAFFT                                    |
| `cast_remap` / `cast_call` | BWA, samtools, bcftools, bgzip, tabix          |


Optional later steps (`cast_extend`, phasing / ploidy) may need additional binaries (e.g. GNU Parallel, exonerate); see the corresponding modules.

---



# Pipeline overview


| Step                       | Input                                   | Output                                                   |
| -------------------------- | --------------------------------------- | -------------------------------------------------------- |
| `cast_assemble`            | `10deduplicated_reads/` + concat probes | `20assemblies/`                                          |
| `cast_retrieve`            | assemblies + split-exon probes          | `30raw_contigs/`, `31exonic_contigs/`                    |
| `cast_analyze`             | `all_hits.tsv`                          | `40aln_orth_par/` (alignments, trees, distances, plots)  |
| `cast_detect`              | hits (+ distances if `-p`)              | `41detected_par/` or `41without_par/`                    |
| `cast_separate`            | customized reference + exonic contigs   | `50pslx/`, `60mafft/`, `70concatenated_exon_alignments/` |
| `cast_remap` → `cast_call` | customized reference + reads            | `100remapped/` (BAMs, VCFs)                              |
| `cast_ploidy`              | remapped BAMs                           | `102ploidy/` (nQuire bins + `lrdmodel.tsv`)              |


Typical local flow:

```bash
python3 ParalogWizard.py cast_assemble  -d DATA -pr probes_concat.fasta -nc 8
python3 ParalogWizard.py cast_retrieve  -d DATA -pe probes_exons.fasta -c -l 75 -s 5 -nc 8
python3 ParalogWizard.py cast_analyze   -d DATA -nc 8
python3 ParalogWizard.py cast_detect    -d DATA -p -mi 4 -ma 10 -nc 8   # or omit -p for no-paralog reference
python3 ParalogWizard.py cast_separate  -d DATA -pc 41detected_par/customized_reference_....fas -i 80 -nc 8
```

---



# Data structure

```
Working_directory/
├── Data_folder/
│   ├── 10deduplicated_reads/
│   ├── 20assemblies/
│   ├── 30raw_contigs/
│   ├── 31exonic_contigs/
│   ├── 40aln_orth_par/
│   ├── 41detected_par/          # with -p
│   ├── 41without_par/           # without -p
│   ├── 50pslx/
│   ├── 60mafft/
│   ├── 70concatenated_exon_alignments/
│   ├── 100remapped/             # remap / call
│   └── 102ploidy/               # nQuire ploidy estimates
├── ParalogWizard/
├── ParalogWizard.py
├── probes_concatenated_exons.fasta
└── probes_separated_exons.fasta
```

Assembly directories under `20assemblies/<sample>/` follow a HybPiper-like per-locus layout for compatibility with existing HybSeq workflows.

---



# Input



### Reads

Trimmed, filtered, deduplicated paired-end reads in `10deduplicated_reads/`, plus `samples_list.txt` (one sample name per line). Names/IDs should use letters, numbers, underscores, and hyphens only.

```
10deduplicated_reads/
├── Genus1-species1_ID.R1.fastq   # or .fastq.gz
├── Genus1-species1_ID.R2.fastq
├── Genus2-species1_ID.R1.fastq
├── Genus2-species1_ID.R2.fastq
└── samples_list.txt
```

`samples_list.txt` example:

```
Genus1-species1_ID
Genus2-species1_ID
Genus2-species2_ID
```



### Probe FASTAs

Two probe files with matching gene names (letters/numbers/hyphens/underscores only):

1. **Concatenated exons** (`-pr` for assemble) — one sequence per gene:

```
>Representative1-Gene1
>Representative1-Gene2
```

1. **Separated exons** (`-pe` for retrieve) — one sequence per exon:

```
>Representative1-Gene1_exon_1
>Representative1-Gene1_exon_2
>Representative1-Gene2_exon_1
```

---



# Local usage

Run from the working directory that contains `ParalogWizard.py`. Use `-h` on any subcommand for the live option list.

### `cast_assemble`

Map reads to concatenated baits (BWA) and assemble (SPAdes). Writes `20assemblies/`.

```bash
python3 ParalogWizard.py cast_assemble -d <data_folder> -pr <concat_probes.fasta> [-nc <cores>]
```


| Option | Description                        |
| ------ | ---------------------------------- |
| `-d`   | Data folder                        |
| `-pr`  | Bait FASTA with concatenated exons |
| `-nc`  | Cores (default 1)                  |




### `cast_retrieve`

Optionally collect SPAdes contigs into `30raw_contigs/` (`-c`), BLAST against split-exon probes, and write exonic contigs plus `all_hits.tsv` / QC under `31exonic_contigs/`.

```bash
python3 ParalogWizard.py cast_retrieve -d <data_folder> -pe <split_exons.fasta> [-c] [-l 75] [-s 5] [-nc <cores>]
```


| Option | Description                                         |
| ------ | --------------------------------------------------- |
| `-d`   | Data folder                                         |
| `-pe`  | Probe FASTA with separated exons                    |
| `-c`   | Collect from `20assemblies/` into `30raw_contigs/`  |
| `-l`   | Minimum probe coverage of BLAST hit (%; default 75) |
| `-s`   | Minimum SPAdes k-mer coverage (default 5)           |
| `-nc`  | Cores (default 1)                                   |




### `cast_analyze`

Per-exon MAFFT alignments and FastTree trees under `40aln_orth_par/`, within-sample pairwise distances, and divergence mixture plots (`pairwise_distances.tsv`, PNGs/SVGs).

```bash
python3 ParalogWizard.py cast_analyze -d <data_folder> [-b <sample> ...] [-nc <cores>]
```


| Option | Description                                                                    |
| ------ | ------------------------------------------------------------------------------ |
| `-d`   | Data folder                                                                    |
| `-b`   | Samples to exclude from divergence estimation (names as in `samples_list.txt`) |
| `-nc`  | Cores (default 1)                                                              |




### `cast_detect`

Build a customized reference.

- **Without** `-p` → `41without_par/` (ortholog-oriented reference, including HybPhyloMaker / HybPiper-oriented FASTA names).
- **With** `-p` → `41detected_par/` (main/para labels from `pairwise_distances.tsv` within `-mi`/`-ma`).

```bash
# No paralog labeling
python3 ParalogWizard.py cast_detect -d <data_folder> [-b <sample> ...] [-nc <cores>]

# With paralog detection (requires cast_analyze)
python3 ParalogWizard.py cast_detect -d <data_folder> -p -mi <min_div> -ma <max_div> [-b <sample> ...] [-nc <cores>]
```


| Option        | Description                                            |
| ------------- | ------------------------------------------------------ |
| `-d`          | Data folder                                            |
| `-b`          | Samples excluded from the written reference            |
| `-p`          | Enable paralog detection                               |
| `-mi` / `-ma` | Min/max % divergence for paralogs (required with `-p`) |
| `-nc`         | Cores (default 1)                                      |




### `cast_separate`

BLAT exonic contigs to the customized reference, then build separated-copy alignments.

```bash
python3 ParalogWizard.py cast_separate -d <data_folder> -pc <customized_reference.fas> -i <min_identity> [-r <sample> ...] [-nc <cores>]
```


| Option | Description                                                                                     |
| ------ | ----------------------------------------------------------------------------------------------- |
| `-d`   | Data folder                                                                                     |
| `-pc`  | Customized reference from `41detected_par/` or `41without_par/`                                 |
| `-i`   | Minimum BLAT identity                                                                           |
| `-r`   | Redlist: pre-WGD / single-copy taxa included in both separated-copy alignments when single-copy |
| `-nc`  | Cores (default 1)                                                                               |




### Downstream (remap / call)

```bash
python3 ParalogWizard.py cast_remap -d <data_folder> -pc <customized_reference.fas> [-e 150] [-nc <cores>]
python3 ParalogWizard.py cast_call  -d <data_folder> [-nc <cores>]
```

---



# Metacentrum usage

```
[server]/
└── home/[logname]/
    ├── Working_directory/
    │   ├── Data_folder/
    │   │   └── 10deduplicated_reads/ ...
    │   ├── ParalogWizard_1a_CastSubmitAssemble.sh
    │   ├── ParalogWizard_1a_CastSubmitAssemble_submitter.sh
    │   ├── ParalogWizard_1b_CastSubmitRetrieve.sh
    │   ├── ParalogWizard_2a_CastSubmitAnalyze.sh
    │   ├── ParalogWizard_2b_CastSubmitDetect.sh
    │   ├── ParalogWizard_3_CastSubmitSeparate.sh
    │   ├── ParalogWizard_4_CastSubmitRemap*.sh / Call / Phase / Ploidy ...
    │   └── ParalogWizard_Settings.cfg
    └── HybSeqSource/
        ├── ParalogWizard/
        ├── ParalogWizard.py
        ├── probes_concatenated_exons.fasta
        ├── probes_separated_exons.fasta
        └── customized_reference_....fas
```



### Settings (`ParalogWizard_Settings.cfg`)


| Setting                                             | Meaning                                                |
| --------------------------------------------------- | ------------------------------------------------------ |
| `data`                                              | Path to data directory under home on the target server |
| `probe_exons_split`                                 | Split-exon probe FASTA (`cast_retrieve -pe`)           |
| `probe_exons_concat`                                | Concatenated probe FASTA (`cast_assemble -pr`)         |
| `exon_length`                                       | Minimum exon length (remap / related steps)            |
| `server`                                            | Metacentrum server name                                |
| `collect_contigs`                                   | `yes` / `no` → `cast_retrieve -c`                      |
| `length_cut`                                        | Probe coverage filter → `-l`                           |
| `spades_cover_cut`                                  | SPAdes k-mer cover filter → `-s`                       |
| `blocklist`                                         | Samples excluded from analyze/detect reference (`-b`)  |
| `redlist`                                           | Pre-WGD / single-copy taxa for separate (`-r`)         |
| `paralogs`                                          | `yes` → detect with `-p`; `no` → without               |
| `paralog_min_divergence` / `paralog_max_divergence` | `-mi` / `-ma`                                          |
| `customized_probes`                                 | Customized reference name for separate (`-pc`)         |
| `minident`                                          | BLAT minimum identity (`-i`)                           |




### Submit scripts


| Script                                                              | Role                                                                             |
| ------------------------------------------------------------------- | -------------------------------------------------------------------------------- |
| `ParalogWizard_1a_CastSubmitAssemble_submitter.sh`                  | Submits one assemble job per sample from `samples_list.txt`                      |
| `ParalogWizard_1a_CastSubmitAssemble.sh`                            | Runs `cast_assemble` (usually via submitter; or `qsub -v sample="[sample]" ...`) |
| `ParalogWizard_1b_CastSubmitRetrieve.sh`                            | Runs `cast_retrieve`                                                             |
| `ParalogWizard_2a_CastSubmitAnalyze.sh`                             | Runs `cast_analyze` (can skip if `paralogs=no`)                                  |
| `ParalogWizard_2b_CastSubmitDetect.sh`                              | Runs `cast_detect`                                                               |
| `ParalogWizard_3_CastSubmitSeparate.sh`                             | Runs `cast_separate` (copy the customized reference into `HybSeqSource/` first)  |
| `ParalogWizard_4_CastSubmitRemap*.sh` / `Call` / `Phase` / `Ploidy` | Downstream variant / ploidy steps                                                |


---



# Logging

Each run writes a single timestamped log, e.g. `ParalogWizard_cast_retrieve_10.Aug.26_12:07.log`. INFO and ERROR always go to that file; `--debug` also writes DEBUG there and prints DEBUG on the console.