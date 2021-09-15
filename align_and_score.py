__author__ = "Josh Myers-Dean"

#from util.util import asterisks_for_residues
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
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
matrix = substitution_matrices.load("BLOSUM62")


def best_align(source: Bio.Seq.Seq, df: pd.DataFrame, records: List[Bio.Seq.Seq], verbose: bool = True) -> List[int]:
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
    return util.get_residues_and_track()

def all_align(source: Bio.Seq.Seq, df: pd.DataFrame, records: List[Bio.Seq.Seq], residues: List[int], seq_id: str, verbose: bool = True, amt: int = 25):
    groups = {
        "A": "1", "I": "1", "L": "1", "V": "1",
        "D": "2", "E":"2", "K": "3", "R": "3", "H": "3",
        "F": "4", "Y": "4", "W": "4", "S": "5", "T": "5",
        "N": "5", "Q": "5", "M": "5", "C": "5", "G": "-1", "P": "-2"
    }
    
    iter = tqdm.tqdm(range(len(records))) if verbose else range(len(records))
    cols = df.columns[3:][::2][:-1]

    all_aligns = {}
    ''' TODO move into function '''
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
            'alignment': align[0],
            'fmt': format_alignment(*align[0])
        } 
    all_aligns = sorted(all_aligns.items(), key=lambda x: x[1]['score'], reverse=True)
    all_aligns = dict(all_aligns)

    count = 0
    for _, val in all_aligns.items():
        if count >= amt:
            break
        src = val['src']
        dest = val['dest']
        fin_inds = util.alignment_index_surgery(residues, src)
        res_score_dict = util.score_based_on_residue(fin_inds, dest, src, groups)
        val.update(res_score_dict)
    all_aligns = sorted(all_aligns.items(), key=lambda v: (v['score'], v['identity'], v['characteristic']))
    all_aligns = dict(all_aligns)
    count = 0
    with open(f'results/{seq_id}.txt', 'w') as f:
        for id, val in all_aligns.items():
            if count > amt:
                break

            length = len(val['src'])
            length = len(val['dest']) if len(val['dest']) > length else length
            alignment_pct = util.alignment_percentage(val['fmt'], length)
            header = f"PDZ domain {id} has an alignment score of {val['score']} ({alignment_pct})%"
            matches = f"\tMatching score is composed of {val['identity']} identity matches and {val['characteristic']} characteristic matches"
            f.write(header + "\n")
            f.write(matches + "\n")
            asterisks_for_residues(residues, f)
            f.write(val['fmt'][: val['fmt'].rfind('c') - 4])
            util.daggers_for_residues(len(val['src']), val['c_inds'], f)
    if verbose:
        with open(f'results/{seq_id}.txt', 'r') as f:
            for line in f.readlines():
                logging.info(line)


if __name__ == '__main__':
    import Bio
    print(Bio.Seq.Seq)
    df = util.clean_csv('data/preds - preds.csv')
    records = util.gen_pdz('data/uniprot_pdz_homosapiens_byDOMAIN_PCLOfix.txt')
    source = SeqIO.parse('fastas/1U3B.fasta', 'fasta')
    src = list(source)[0].seq
    #residues = best_align(src, df, records)
    residues = [109, 111, 116, 119, 148, 152, 156]
    residues = [x-1 for x in residues]
    all_align(src, df, records, residues, '1U3B', True)
    
