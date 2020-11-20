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
  * [MAFFT 6.9 or later](https://mafft.cbrc.jp/alignment/software/)
___
**Input and prior requirements**


**Usage**

ParalogWizard_CastCollect.py "folder with HybPiper results ('data_HybPiper')" [folder for data output as in HybPhyloMaker ('data')]

ParalogWizard_CastConvert.py [folder for data output as in HybPhyloMaker ('data')] [probe file with separated exons ('probe_HP_exons_split')] [threshold for length cover of BLAST hits ('length_cut')] "threshold for k-mer cover of contigs assembled by SPAdes ('spades_cover_cut')" 



