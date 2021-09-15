__author__ = "Josh Myers-Dean"

from Bio import SeqIO
from Bio import pairwise2
import Bio
from Bio.pairwise2 import format_alignment
import pandas as pd
import os
from Bio.Align import substitution_matrices
import tqdm
from util import util
from typing import List, NoReturn
import logging
logging.basicConfig(level=logging.INFO, format='%(message)s')
matrix = substitution_matrices.load("BLOSUM62")


def best_align(source: Bio.Seq.Seq, df: pd.DataFrame, records: List[Bio.Seq.Seq], verbose: bool = True) -> List[int]:
    '''
    Find best aligment for input PDZ domain out of 140 labeled domains.

    params:
        source: Input PDZ domain
        df: dataframe of labelled data
        records: list of human pdz domains
        verbose: have verbose terminal output 
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

def all_align(source: Bio.Seq.Seq, df: pd.DataFrame, records: List[Bio.Seq.Seq], residues: List[int], seq_id: str, verbose: bool = True, amt: int = 25) -> NoReturn:
    '''
    aligns input PDZ domain with all 272 human PDZ domains. Computes alignment score, identity, and characteristic matches.

    params
        source: Input PDZ domain
        df: dataframe of labelled data
        records: list of human pdz domains
        residues: list of conserved residue positions in input domain
        seq_id: name of the sequence, used for saving the file
        verbose: have verbose terminal output
        amt: amount of matches to show
    '''
    groups = {
        "A": "1", "I": "1", "L": "1", "V": "1",
        "D": "2", "E":"2", "K": "3", "R": "3", "H": "3",
        "F": "4", "Y": "4", "W": "4", "S": "5", "T": "5",
        "N": "5", "Q": "5", "M": "5", "C": "5", "G": "-1", "P": "-2"
    }
    
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
            'alignment': align[0],
            'fmt': format_alignment(*align[0])
        } 
    all_aligns = sorted(all_aligns.items(), key=lambda x: x[1]['score'], reverse=True)
    all_aligns = dict(all_aligns)

    for _, val in all_aligns.items():
        src = val['src']
        dest = val['dest']
        fin_inds = util.alignment_index_surgery(residues, src)
        res_score_dict = util.score_based_on_residue(fin_inds, dest, src, groups)
        val.update(res_score_dict)

    all_aligns = sorted(all_aligns.items(), key=lambda v: (v[1]['score'], v[1]['identity'], v[1]['characteristic']), reverse=True)
    all_aligns = dict(all_aligns)
    count = 0
    os.makedirs('results', exist_ok=True)
    with open(f'results/{seq_id}.txt', 'w') as f:
        for id, val in all_aligns.items():
            if count >= amt:
                break
            length = len(val['src'])
            alignment_pct = util.alignment_percentage(val['fmt'], length) * 100
            header = f"PDZ domain {id} has an alignment score of {val['score']} ({alignment_pct:.2f})%"
            matches = f"\tMatching score is composed of {val['identity']} identity matches and {val['characteristic']} characteristic matches"
            f.write(header + "\n")
            f.write(matches + "\n")
            util.asterisks_for_residues(residues, f)
            f.write(val['fmt'][: val['fmt'].rfind('c') - 4])
            util.daggers_for_residues(len(val['src']), val['c_inds'], f)
            count = count + 1
    if verbose:
        with open(f'results/{seq_id}.txt', 'r') as f:
            logging.info(f.read())


if __name__ == '__main__':
    args = util.get_args()
    df = util.clean_csv('data/preds - preds.csv')
    records = util.gen_pdz('data/uniprot_pdz_homosapiens_byDOMAIN_PCLOfix.txt')
    source = SeqIO.parse(args.f, 'fasta')
    fname = args.f.split('/')[-1].replace('.fasta', '')
    src = list(source)[0].seq

    if args.p is not None:
        residues = util.convert(args.p)
    else:
        residues = best_align(src, df, records)
    all_align(src, df, records, residues, fname, args.v, args.m)
    
