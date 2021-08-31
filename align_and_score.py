__author__ = "Josh Myers-Dean"

from Bio import SeqIO
from Bio import pairwise2
import Bio
from Bio.pairwise2 import format_alignment
import pandas as pd
import os
import argparse
import sys
import numpy as np
from Bio.Align import substitution_matrices
import tqdm
from util import util
from typing import List
matrix = substitution_matrices.load("BLOSUM62")


def best_align(source: Bio.Seq.Seq, df: pd.DataFrame, records: List[Bio.Seq.Seq], verbose: bool = True):
    '''
    TODO
    '''
    records = [i for i in records if i.id in df['pdz_domain'].to_list()]
    inds = [df[x] for x in [str(i) for i in range(1,8)]]
    iter = tqdm.tqdm(range(len(records))) if verbose else range(len(records))
    all_aligns = {}

    for s in iter:
        
        id_ = records[s].id
        seq = records[s].seq
        align = pairwise2.align.globalds(source, seq, matrix, -11, -1)
        src_aligned = align[0].seqA
        dest_aligned = align[0].seqB
        all_aligns[id_] = {
            'src': src_aligned,
            'dest': dest_aligned,
            'score': align[0].score,
            'alignment': align[0]
        }
    all_aligns = sorted(all_aligns.items(), key=lambda x: x[1]['score'], reverse=True)
    all_aligns = dict(all_aligns)
    util.pretty_print_first_match(all_aligns, df, inds, list(all_aligns.keys())[0])
    return util.get_residues_and_start()

def all_align(source: Bio.Seq.Seq, df: pd.DataFrame, records: List[Bio.Seq.Seq], verbose: bool = True, amt: int = 25):
    pass



if __name__ == '__main__':
    df = util.clean_csv('data/preds - preds.csv')
    records = util.gen_pdz('data/uniprot_pdz_homosapiens_byDOMAIN_PCLOfix.txt')
    source = SeqIO.parse('fastas/1U3B.fasta', 'fasta')
    src = list(source)[0].seq
    print(best_align(src, df, records))

    
