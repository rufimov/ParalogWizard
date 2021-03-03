<pre> 
  _____                _          __          ___                  _ 
 |  __ \              | |         \ \        / (_)                | |              /\
 | |__) |_ _ _ __ __ _| | ___   __ \ \  /\  / / _ ______ _ _ __ __| |             /  \
 |  ___/ _` | '__/ _` | |/ _ \ / _` \ \/  \/ / | |_  / _` | '__/ _` |            |    |
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
                                                     '-   |_         -</pre>
---
# Dependencies
  * [Python 3.6 or later](https://www.python.org/downloads/)
  * [BIOPYTHON 1.77 or later](https://biopython.org/wiki/Download)
  * [BLAST command line tools 2.2.30 or later](https://www.ncbi.nlm.nih.gov/books/NBK131777/#_Blast_ReleaseNotes_BLAST_2_2_30_October_)
  * [NumPy](https://numpy.org)
  * [Scikit-learn](https://scikit-learn.org/stable/user_guide.html)
  * [SciPy](https://www.scipy.org)
  * [Matplotlib](https://matplotlib.org)
  * [MAFFT 6.9 or later](https://mafft.cbrc.jp/alignment/software/)
  * [FastTree 2.0.0 or later](http://www.microbesonline.org/fasttree/)
  * [BLAT](http://genome.ucsc.edu/goldenPath/help/blatSpec.html)
  * [SPAdes](http://cab.spbu.ru/software/spades/)
  * [GNU Parallel](http://www.gnu.org/software/parallel/) 
  * [BWA](https://github.com/lh3/bwa)
  * [samtools 1.2 or later](https://github.com/samtools/samtools)

---
# Data structure

<pre>
Working_directory
    |_____Data_folder
    |         |_____10deduplicated_reads
    |         |_____20assemblies
    |         |_____30raw_contigs
    |         |_____31exonic_contigs
    |         |_____40aln_orth_par
    |         |_____41detected_par
    |         |_____50pslx
    |         |_____60mafft
    |         |_____70concatenated_exon_alignments
    |_____ParalogWizard
    |_____ParalogWizard.py
    |_____probes_concatenated_exons.fasta
    |_____probes_separated_exons.fasta
</pre>
    
---
# Input

* Trimmed, filtered and deduplicated pair-end reads are store in 10deduplicated_reads folder together with list of samples called ```samples_list.txt```. All names and IDs shoud contain only letters and numbers. 

<pre>
10deduplicated_reads
    |_____Genus1-species1_ID.R1.fastq
    |_____Genus1-species1_ID.R2.fastq
    |_____Genus2-species1_ID.R1.fastq
    |_____Genus2-species1_ID.R2.fastq
    |_____Genus2-species2_ID.R1.fastq
    |_____Genus2-species2_ID.R2.fastq
    ...
    |_____samples_list.txt
</pre>

List of samples must have names of the samples which correspond to fastq files, one per each line with empy line in the end.
```    
    Genus1-species1_ID
    Genus2-species1_ID
    Genus2-species2_ID
    ...
    
```

* Two probe files: one containing sequences which are concatendated exons, another - with the same sequences separated to exons. Both files should have same names for the same gene. All names shoud contain only letters and numbers.
```
    >Representative1-Gene1
    >Representative1-Gene2
    >Representative2-Gene3
    ...
 ```

```
    >Representative1-Gene1_exon_1
    >Representative1-Gene1_exon_2
    >Representative1-Gene2_exon_1
    >Representative2-Gene3_exon_2
    >Representative2-Gene3_exon_4
    ...
```

---
# Local usage


```python3 ParalogWizard.py cast_retrieve -d <folder with data> -pe <probe file with separated exons> -l <threshold for length cover of BLAST hits> -s <threshold for k-mer cover of contigs assembled by SPAdes> [-nc <number of cores>] [-c]```

Collects contigs assembled by SPAdes to folder ```30raw_contigs```. Matches retrieved contigs to the probe file with individual separated exons with BLAST, extracts exonic contigs according to hits and stores them in folder ```31exonic_contigs``` within the main folder with ParalogWizards results. Hit tables and statistics with asumed number of copies for each locus per sample are saved alongside.

```python3 ParalogWizard.py cast_analyze -d <folder with data> [-b <list of taxa excluded from paralog divergence estimation>] [-nc <number of cores>]```

```python3 ParalogWizard.py cast_detect -d <folder with data> -pe <probe file with separated exons> [-b <list of taxa excluded from new reference creation>] [-p -mi <minimum paralog divergence> -ma <maximum paralog divergence>] ``` 

```python3 ParalogWizard.py cast_separate -d <folder with data> -pp <probe file with separated paralogs> -i <minimum identity for BLAT> [-r <list of taxa excluded from paralogs separation, ie included in all alignments in case of >]```


---
# Metacentrum usage




