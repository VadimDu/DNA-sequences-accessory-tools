
# title: random_sampling_fasta.py
# description: Tool for random sampling (without replacment) of sequences in a fasta file
# author: Vadim (Dani) Dubinsky [dani.dubinsky@gmail.com]
# date: 03/05/2021
# required modules: Biopython

import argparse
import os
import sys
from collections import OrderedDict
import random

try:
    from Bio import SeqIO
except:
    sys.exit("This program requires Python3 Bio module, please install it")

def main():

    args = args_setup()

    path_exists(args.fasta)

    rand_sampling_fasta(args.fasta, args.output, args.sequence_num)


def rand_sampling_fasta(in_fasta, out_fasta, sequence_num):
    ''' Random sampling (without replacment) of sequences in a fasta file according to a specified number of sequences to sample
        Parameters
        ----------
        in_fasta : input fasta file (path)
        out_fasta : output file for the sampled subset fasta file (path)
        sequence_num : required number of sequences to sample from the input fasta 
        Returns
        -------
        None.
        '''

    frag_count = 0
    frags_odict = OrderedDict()

    for frag in SeqIO.parse(in_fasta, "fasta"):
        frag_count += 1
        frags_odict[str(frag.description)] = frag.seq

    #Random sampling of x_rand_frags fragments, without replacment:
    if (sequence_num < frag_count):
        frags_odict_rand = OrderedDict(random.sample(frags_odict.items(), sequence_num))
    else:
        sys.exit("Please specify a smaller number to sample, as your input fasta contains only: " + str(frag_count) + " sequences")

    with open(out_fasta, 'w') as f:
        for keys,vals in frags_odict_rand.items():
            f.write(">" + keys + "\n" + str(vals) + "\n")


def path_exists(file):
    if (not os.path.exists(file)):
        sys.exit("File path error: " + file + " does not exist!")


def args_setup():
    '''Command line arguments parsing function'''
    parser = argparse.ArgumentParser(description="Tool for random sampling (without replacment) of fasta files")
    parser.add_argument("fasta", help="fasta file with nucleotide/aa sequences")
    parser.add_argument("-s", "--sequence_num", metavar="<sequence number (int)>", type=int, required=True,
                        help="the required number of sequences to sample from the input fasta")
    parser.add_argument("-o", "--output", metavar='<path-to-output>', required=True,
                        help="path to output file name for the sampled subset of the original fasta")
    return parser.parse_args()


if __name__ == "__main__":
    main()
