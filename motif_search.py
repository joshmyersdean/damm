__author__ = 'Josh Myers-Dean'

from Bio import SeqIO
from util import config
import argparse
import logging
from typing import NoReturn
import os
logging.basicConfig(level=logging.INFO, format='%(message)s')


def get_args() -> config.CfgNode:
    '''
    Grabs and parses config file
    '''
    parser = argparse.ArgumentParser(description='Motif Matcher Search')
    parser.add_argument(
        '--config',
        type=str,
        default='config/example.yaml',
        help='config file')
    args = parser.parse_args()
    assert args.config is not None
    cfg = config.load_cfg_from_cfg_file(args.config)
    return cfg


def position_mapping(pos: str) -> int:
    '''
    Helper function to motif positions to program positions.
    args
            pos: position of motif residue
    '''
    positions = {
        'p0': -1,
        'p-1': -2,
        'p-2': -3,
        'p-3': -4,
        'p-4': -5,
        'p-5': -6,
    }
    return positions[pos]


def check_verbose(args) -> bool:
    try:
        return args.verbose
    except KeyError:
        return False


def search(args: config.CfgNode) -> NoReturn:
    src = args.motif
    max_subs = args.max_subs
    verbose = check_verbose(args)
    os.makedirs(args.save_dir, exist_ok=True)
    with open(f"{args.save_dir}/{args.file_name}", 'w') as f:
        f.write("subs \t \t\tmatch \t\t matchSeq \t src\n")
        f.write('-' * 65)
        f.write('\n')
        with open(args.proteome) as handle:
            for record in SeqIO.parse(handle, "fasta"):
                free_positions = [-6, -5, -4, -3, -2, -1]
                constrained_pos = []
                constraints_met = True
                subs = 0
                for arg, val in args.items():
                    if arg.startswith('p') and arg != 'proteome':
                        position = position_mapping(arg)
                        constrained_pos.append(position)
                        if record.seq[position] in val:
                            continue
                        else:
                            constraints_met = False
                if not constraints_met:
                    continue
                free_positions = [
                    i for i in free_positions if i not in constrained_pos]
                for pos in free_positions:
                    if record.seq[pos] != src[pos]:
                        subs = subs + 1
                        if subs > max_subs:
                            constraints_met = False
                if constraints_met:
                    f.write(
                        f"{subs} \t {record.id} \t {record.seq[-6:]} \t {src}\n")
    if verbose:
        with open(f"{args.save_dir}/{args.file_name}", 'r') as f:
            logging.info(f.read())


if __name__ == '__main__':
    search(get_args())
