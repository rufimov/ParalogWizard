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



===

___
**Dependencies**
  * [BLAST command line tools 2.2.30 or later](https://www.ncbi.nlm.nih.gov/books/NBK131777/#_Blast_ReleaseNotes_BLAST_2_2_30_October_)
  * [Python 3.6 or later](https://www.python.org/downloads/)
  * [BIOPYTHON 1.77 or later](https://biopython.org/wiki/Download)
  * [NumPy](https://numpy.org)
  * [Scikit-learn](https://scikit-learn.org/stable/user_guide.html)
  * [SciPy](https://www.scipy.org)
  * [Matplotlib](https://matplotlib.org)
  * [MAFFT 6.9 or later](https://mafft.cbrc.jp/alignment/software/)
  * [FastTree 2.0.0 or later](http://www.microbesonline.org/fasttree/)
___
**Input and prior requirements**

Genus-species_ID

Representative-Gene_exon_num




**Local usage**

```python3 ParalogWizard.py cast_collect -c <folder with HybPiper results> -d <folder for data output as in HybPhyloMaker>```

Collects contigs assembled by SPAdes within HybPiper and stores in folder 'HybPiper_contigs' within the main folder with ParalogWizards results.


```python3 ParalogWizard.py cast_retrieve -d <folder for data output as in HybPhyloMaker> -pe <probe file with separated exons> -l <threshold for length cover of BLAST hits> -s <threshold for k-mer cover of contigs assembled by SPAdes> [-nc <number of cores>]```

Matches retrieved contigs to the probe file with individual separated exons with BLAST, extracts exonic contigs according to hits and stores them in folder 'exons/40contigs' within the main folder with ParalogWizards results. Hit tables and statistics with asumed number of copies for each locus per sample are saved alongside.

```python3 ParalogWizard.py cast_analyze -d <folder for data output as in HybPhyloMaker> [-b <list of taxa excluded from paralog divergence estimation>] [-nc <number of cores>]```

```python3 ParalogWizard.py cast_create -d <folder for data output as in HybPhyloMaker> -pe <probe file with separated exons> [-b <list of taxa excluded from new reference creation>] [-p -mi <minimum paralog divergence> -ma <maximum paralog divergence>] ``` 

```python3 ParalogWizard.py cast_correct -d <folder for data output as in HybPhyloMaker> -pp <probe file with separated paralogs> -i <minimum identity for BLAT> [-r <list of taxa excluded from paralogs separation, ie included in all alignments in case of >]```



**Metacentrum usage**



