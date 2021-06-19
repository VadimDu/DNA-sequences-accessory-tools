
# title: handle_none_ATGC_fasta.py
# description: Tool for handling none-ATGC nucleotides in a fasta file of nucleotide sequences
# author: Vadim (Dani) Dubinsky [dani.dubinsky@gmail.com]
# date: 05/05/2021
# required modules: Biopython

import argparse
import os
import sys
from collections import OrderedDict
import random
import re

try:
    from Bio import SeqIO
except:
    sys.exit("This program requires Python3 Bio module, please install it")


def main():
    
    args = args_setup()
    
    path_exists(args.fasta)
    
    handle_none_ATGC(args.fasta, args.output)


def handle_none_ATGC(in_fasta, out_fasta):
    ''' Handle none-ATGC nucleotides in a fasta file by replacing them with a random nucleotides
    Parameters
    ----------
    in_fasta : input fasta file (path)
    out_fasta : output file for the corrected non-ATGC fasta (path)
    Returns
    -------
    None.
    '''
    seq_odict = OrderedDict()
    ambig_count = 0
    ATGC_bases = ["A", "T", "G", "C"]

    for sequence in SeqIO.parse(in_fasta, "fasta"):
        # check if an ambiguous base is found
        if re.search("[^ATGC]", str(sequence.seq)):
            new_seq = re.sub("[^ATCG]", random.choice(ATGC_bases), str(sequence.seq))
            assert len(sequence.seq) == len(new_seq)
            seq_odict[str(sequence.description)] = new_seq
        else:
            seq_odict[str(sequence.description)] = sequence.seq
    
        m = re.search("[^ATGC]", str(sequence.seq))
        if m:
            ambig_count += 1
            print("ambiguous base \"%s\" found!" % (m.group()))
    
    print("Total sequences with ambiguous bases found:", ambig_count)
    
    with open(out_fasta, 'w') as f:
        for keys, vals in seq_odict.items():
            f.write(">" + keys + "\n" + str(vals) + "\n")    
    
    
def path_exists(file):
        if (not os.path.exists(file)):
            sys.exit("File path error: " + file + " does not exist!")


def args_setup():
    '''Command line arguments parsing function'''
    parser = argparse.ArgumentParser(description="Tool for handling none-ATGC nucleotides in a fasta file by replacing them with a random nucleotides (By Vadim-Dani Dubinsky dani.dubinsky@gmail.com)")
    parser.add_argument("fasta", help = "fasta file with nucleotide contigs (also accepts amino-acid sequences)")
    parser.add_argument("-o", "--output", metavar = '<path-to-output>', required = True, help = "path to output file name for the corrected non-ATGC fasta")
    return (parser.parse_args())


if __name__ == "__main__":
    main()
