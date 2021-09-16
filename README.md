# DAMM - Domain Analysis and Motif Matcher


## Requirements
Python 3.6 or greater
```bash
pip install -r requirements.txt
```

## Domain Analysis
Domain analysis aims to find similar PDZ domains to an input domain based off of total alignment, identity matches at conserved residues, and characteristic matches at conserved residues. We supply a list of 140 labelled human PDZ domains to help locate conserved residues in the input domain.

If you do not know the located of the conserved residues of your input domain, run the program as such:
```bash
python align_and_score.py -f {FASTA FILE OF YOUR DOMAIN}
```
to be presented with the top alignment, along with the conserved residues of that alignment. You will then be prompted to input the positions of your 7 residues, followed by the starting position of your sequence (e.g. if your sequence starts at 715 and your first residue is at 733, use those numbers).  

If you wish to bypass this initial matching step you can simply provide the **0-indexed** positions of your sequence as such:

```bash
$ python align_and_score.py -f {FASTA} --p {p1,p2,p3,p4,p5,p6,p7}

nput aligned with APBA1-2, alignment score: 258.0 (41%)
         10        20        30        40        50        60        70        80        90        100       110       120       130       140       150       160       170       180  
EFKDVFIEKQKGEILGVVIVESGWGSILPTVIIANMMHGGPAEKSGKLNIGDQIMSINGTSLVGLPLSTCQSIIKGLKNQSRVKLNIVRCPPVTTVLIRRPDLRYQLGFSVQNGIICSLMRGGIAERGGVRVGHRIIEINGQSVVATPHEKIVHILSNAVGEIHMKTMPAAMYRLLTAQEQPVYI
                                                                                              ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||               
----------------------------------------------------------------------------------------------TVLIRRPDLRYQLGFSVQNGIICSLMRGGIAERGGVRVGHRIIEINGQSVVATPHEKIVHILSNAVGEIHMKTMPA---------------
                                                                                                             * *    *  *                            *   *   *                            

Input 7 Positions (1-indexed):
        ex: 1,2,3,4,5,6,7
Input positions: 126,128,133,136,151,165,169
Input starting residue position (1 indexed):
Starting position: 17
```

We also offer the following flags
- `--v`: show verbose terminal output
- `--m`: number of matches to show

### Results
Results are saved by default in the `results` folder.

### Example Usage
```bash
$ python align_and_score.py -f fastas/1U3B.fasta --v True --m 50 --p 109,111,116,119,148,152,156
```

### Retreiving Fasta Files
We supply a shell script to download fasta files from the PDB. Simply call the script with the PDB ID of your desired sequence.
#### Example Usgae
```bash
$ bash get_fasta.sh 1U3B
```

## Motif Matcher
To use our implementation of Motif Matcher, simply create or a modify a config file. We provide an example in the `configs` directory. Note that you must specify a proteome file and motif.

### Example Usage
```bash
$ python motif_search.py --config {config name}

subs                    match            matchSeq        src
-----------------------------------------------------------------
2        tr|A9V6G5|A9V6G5_MONBE          HRESTV          HSETAL
2        tr|A9V7Z4|A9V7Z4_MONBE          EDETAL          HSETAL
2        tr|A9UWK2|A9UWK2_MONBE          LSESRI          HSETAL
2        tr|A9V752|A9V752_MONBE          TSESRL          HSETAL
2        tr|A9UXE1|A9UXE1_MONBE          QDETAL          HSETAL
1        tr|A9UNL4|A9UNL4_MONBE          HSETTF          HSETAL
2        tr|A9V727|A9V727_MONBE          SSESRL          HSETAL
2        tr|A9UP37|A9UP37_MONBE          LSGTAI          HSETAL
2        tr|A9UP44|A9UP44_MONBE          QSESRL          HSETAL
2        tr|A9VBW1|A9VBW1_MONBE          FSHTAL          HSETAL
2        tr|A9V3B6|A9V3B6_MONBE          RSLTAI          HSETAL
2        tr|A9UWU2|A9UWU2_MONBE          YSETYV          HSETAL
2        tr|A9VAF6|A9VAF6_MONBE          LSESVV          HSETAL
2        tr|A9VBH4|A9VBH4_MONBE          TSESVL          HSETAL
2        tr|A9UWP5|A9UWP5_MONBE          GSESSV          HSETAL
2        tr|A9V5L4|A9V5L4_MONBE          HSKSFL          HSETAL
2        tr|A9UTX9|A9UTX9_MONBE          ESESMV          HSETAL
2        tr|A9V4H3|A9V4H3_MONBE          ISESCL          HSETAL
2        tr|A9VBG8|A9VBG8_MONBE          VSTSAL          HSETAL
2        tr|A9UWL5|A9UWL5_MONBE          EEESAV          HSETAL
```

## Development
All development was done on MacOS and Ubuntu 20.04 and is not tested on Windows.