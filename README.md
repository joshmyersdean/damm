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
python align_and_score.py -f {FASTA} --p {p1,p2,p3,p4,p5,p6,p7}
```

We also offer the following flags
- `--v`: show verbose terminal output
- `--m`: number of matches to show

### Example Usage
```bash
python align_and_score.py -f fastas/1U3B.fasta --v True --m 50 --p 109,111,116,119,148,152,156
```

### Retreiving Fasta Files
We supply a shell script to download fasta files from the PDB. Simply call the script with the PDB ID of your desired sequence.
#### Example Usgae
```bash
bash get_fasta.sh 1U3B
```

## Development
All development was done on MacOS and Ubuntu 20.04 and is not tested on Windows.

## Motif Matcher
TODO