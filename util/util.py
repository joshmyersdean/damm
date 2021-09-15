__author__ = 'Josh Myers-Dean'

import pandas as pd
from typing import List, Dict, Any, NoReturn, Tuple
import Bio
import argparse
from Bio.pairwise2 import format_alignment
import io

def convert(inp: str) -> List[int]:
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
    if not inp: 
        return []
    pages = []
    comma_separated = []
    comma_separated = inp.split(",") 
    for item in comma_separated:
        if "-" in item:
            a = item.split("-")
            pages.extend(range(int(a[0]),int(a[1])+1))
        else:
            pages.append(int(item))
    return pages


def clean_csv(fname: str) -> pd.DataFrame:
    '''
    Cleans up columns in dataframe for downstream use.

    params
        fname: name of csv file
    '''

    df = pd.read_csv(fname)
    df = df.drop(columns=df.columns[-1])
    df['pdz_domain'] = df['pdz_domain'].apply(lambda x: x.replace('>', ''))
    df = df.dropna().reset_index(drop=True)
    return df

def gen_pdz(fname: str) -> List[Bio.Seq.Seq]:
    '''
    generates list of sequences of human PDZ domains.
    params
        fname: name of uniprot file.
    '''

    records = []
    with open(fname, 'r') as handle:
        for record in Bio.SeqIO.parse(handle, 'fasta'):
            records.append(record)
    return records

def pretty_print_first_match(alignments: Dict[str, Dict[str, Any]], df: pd.DataFrame, inds: List[pd.Series], id: str) -> NoReturn:
    '''
    Syntatic sugar helper method for showing user top scoring match.
    params
        alignments: Dictionary of alignments with score, sequences, alignment.
        df: pandas dataframe of labelled domains
        inds: labelled indices of pdz domains
        id: name of top scoring PDZ domain
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

def get_residues_and_track() -> List[int]:
    '''
    Gets user input of conserved residues and shifts to 0 indexed values.
    returns
        list of residues
    raises
        AssertionError: must enter 7 positions, must be sorted, start must be positive
    '''
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
        print("Input Error.")
        return get_residues_and_track()
    positions = list(map(lambda x: x-start, positions))
    return positions

def alignment_index_surgery(weights: List[int], src: Bio.Seq.Seq) -> List[int]:
    '''
    Get positions of conserved residues after alignment to account for gaps.
    params
        weights: List of conserved residues
        src: input PDZ domain
    returns
        List of aligned residues
    '''
    fin_inds = []

    for i in weights:
        tmp_inds = []
        track = 0
        dont_break = True
        for j in range(len(src)):
            if j < track:
                continue
            while src[track] == '-' and dont_break:
                if track == len(src)-1:
                    dont_break = False
                    break
                track += 1
            tmp_inds.append(track)
            track += 1
        fin_inds.append(tmp_inds[i])
    return fin_inds


def score_based_on_residue(fin_inds: List[int], dest: Bio.Seq.Seq, src: Bio.Seq.Seq, groups: Dict[str, int]) -> Dict[str, int]:
    '''
    Calculate number of identity and characteristic matches, as well as characteristic match positions.
    params
        find_inds: aligned indices for residues
        dest: destination PDZ domain
        src: source PDZ domain
        groups: Characteristic residue groups
    returns
        dictionary of number of identity matches, characteristic matches, and characteristic indices
    '''
    identity = 0
    characteristic = 0
    characteristic_inds = []
    for ind in fin_inds:
        if ind > len(dest):
            continue
        if src[ind] == '-' or dest[ind] == '-':
            continue
        elif src[ind] == dest[ind]:
            identity += 1
        elif groups[src[ind]] == groups[dest[ind]]:
            characteristic += 1
            characteristic_inds.append(ind)

    ret_dict = {
        'identity': identity,
        'characteristic': characteristic,
        'c_inds': characteristic_inds
    }  
    return ret_dict

def alignment_percentage(alignment: str, length: int) -> float:
    '''
    Calculates percentage of alignment matches/
    params
        alignment: formatted alignment
        length: length of alignment including gaps
    returns
        alignment percentage
    '''
    numerator = alignment.count('|')
    return numerator / length

def asterisks_for_residues(residues: List[int], fobj: io.TextIOWrapper) -> NoReturn:
    '''
    Writes asterisks where conserved residues exist to file.
    params
        residues: list of conserved residues
        fobj: file to write to
    '''
    track = 0
    for i in residues:
        while (track < i):
            fobj.write(" ")
            track = track + 1
        fobj.write('*')
        track = track + 1
    fobj.write("\n")

def daggers_for_residues(len_seq: int, char_inds: List[int], fobj: io.TextIOWrapper):
    '''
    Writes daggers where characteristic matches exit to file.
    params
        residues: list of conserved residues
        char_inds: list of indices where characteristic matches occur
        fobj: file to write to 
    '''
    dagger = u"\u2020"
    fobj.write("\n")
    for idx in range(len_seq):
        if idx in char_inds:
            fobj.write(dagger)
        else:
            fobj.write(" ")
    fobj.write('\n\n')

def get_args() -> argparse.Namespace:
    '''
    Retreives command line arguments.
    returns
        namespace with arguments
    '''
    parser = argparse.ArgumentParser(description='Match a given FASTA sequence with one of the 272 Human PDZ Domains',
                                    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', type=str, help='Fasta file containing a single sequence for matching')
    parser.add_argument('--m', type=int, default=25, help='amount  of matches to show')
    parser.add_argument('--p', type=str, default=None, help='0-indexed conserved residue positions')
    parser.add_argument('--v', type=bool, default=True, help='Verbose output: Display file to terminal')
    return parser.parse_args()