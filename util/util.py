__author__ = 'Josh Myers-Dean'

import pandas as pd
from typing import List, Dict, Any, NoReturn, Tuple
import Bio
from Bio.pairwise2 import format_alignment

def convert(inp):
    """
    https://gist.github.com/amjith/673824
    * Get the input from user.
    * Parse the input to extract numbers
    ** Split by comma
    *** Each item in the list will then be split by '-'
    **** Populate the number between a-b using range(a,b)
    >>> convert("")
    []
    >>> convert("1")
    [1]
    >>> convert("1,2")
    [1, 2]
    >>> convert("1,2-5")
    [1, 2, 3, 4, 5]
    >>> convert("1-3,2-5,8,10,15-20")
    [1, 2, 3, 2, 3, 4, 5, 8, 10, 15, 16, 17, 18, 19, 20]
    """
    if not inp:  # empty string must return an empty list
        return []
    pages = []
    comma_separated = []
    comma_separated = inp.split(",")     # split the input string based on comma delimitation
    for item in comma_separated:
        if "-" in item:
            a = item.split("-")
            pages.extend(range(int(a[0]),int(a[1])+1))
        else:
            pages.append(int(item))
    return pages


def clean_csv(fname: str) -> pd.DataFrame:
    '''
    TODO
    '''
    df = pd.read_csv(fname)
    df = df.drop(columns=df.columns[-1])
    df['pdz_domain'] = df['pdz_domain'].apply(lambda x: x.replace('>', ''))
    df = df.dropna().reset_index(drop=True)
    return df

def gen_pdz(fname: str) -> List[Bio.Seq.Seq]:
    '''
    TODO
    '''
    records = []
    with open(fname, 'r') as handle:
        for record in Bio.SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

def pretty_print_first_match(alignments: Dict[str, Dict[str, Any]], df: pd.DataFrame, inds: List[pd.Series], id: str) -> NoReturn:
    '''
    TODO
    '''
    best_match = alignments[id]
    alignment = format_alignment(*best_match['alignment'])
    ticks = [" " for i in range(len(best_match['dest']))]
    index = df[df['pdz_domain'] == id].index.values[0]
    print(f"Input aligned with {id}, alignment score: {best_match['score']} ({int(100*alignment.count('|') / len(best_match['dest']))}%)")
    for i in range(0,len(ticks)):
        if i == 0:
            print(" ", end="")
            continue
        # if i == len(ticks)-1:
        #     print("", end="\n")
        if i >= 100:
            if (i-1) % 10 == 0 or (i-2) % 10 == 0:
                continue
        else:
            if (i-1) % 10 == 0:
                continue
        if i % 10 == 0:
            if i < 100:
                print(f"{i // 1}", end="")
            else:
                print(f"{i // 1}", end="")
        else:
            print(" ", end="")
    
    for p in range(7):
        absolute_pos = inds[p][index] - df['1st residue'][index]
        ticks[absolute_pos] = '*'
    print("\n"+alignment[:alignment.rfind("c")-4], end='\n')
    print(''.join(ticks))
    print()

def get_residues_and_start() -> List[int]:
    print("Input 7 Positions (1-indexed):\n\tex: 1,2,3,4,5,6,7")
    positions = input("Input positions: ").strip()
    positions = convert(positions)
    print("Input starting residue position (1 indexed):")
    start = input("Starting position: ").strip()
    start = int(start)
    try:
        assert len(positions) == 7, 'Must enter 7 positions.'
        assert sorted(positions) == positions, 'Positions must be sorted'
        assert start > 0, 'Starting index must be positive'
    except AssertionError:
        return get_residues_and_start()
    positions = list(map(lambda x: x-start, positions))
    return positions

    