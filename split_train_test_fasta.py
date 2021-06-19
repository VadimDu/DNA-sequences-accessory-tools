
# title: random_sampling_fasta.py
# description: Tool for splitting a fasta files into 2 separates files - train set & test set, based on a test size proportion
# author: Vadim (Dani) Dubinsky [dani.dubinsky@gmail.com]
# date: 05/05/2021
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
    
    split_train_test(args.fasta, args.output_train, args.output_test, args.size_testset)


def split_train_test(in_fasta, out_fasta_train, out_fasta_test, test_size):
    ''' Split a fasta files into 2 separates files - train set & test set, based on a test size proportion
        Parameters
        ----------
        in_fasta : input fasta file (path)
        out_fasta_train : output file for train set (path)
        out_fasta_test : output file for test set (path)
        test_size : required size (proportion) of the test set
        Returns
        -------
        None.
        '''

    seq_odict = OrderedDict()
    seq_count = 0
    
    for sequence in SeqIO.parse(in_fasta, "fasta"):
        seq_count += 1
        seq_odict[str(sequence.description)] = sequence.seq
    
    #the actual number of test set sequences
    test_key_count = int((len(seq_odict.keys()))*test_size)

    #random sampling of test_key_count sequences (which is test_size % of the whole data), without replacment:
    test_keys = random.sample(seq_odict.keys(), test_key_count)
    train_keys = [s for s in seq_odict.keys() if s not in test_keys]

    #generate the splitted train-test dict:
    test_odict = OrderedDict((key, seq_odict[key]) for key in test_keys if key in seq_odict)
    train_odict = OrderedDict((key, seq_odict[key]) for key in train_keys if key in seq_odict)

    print ("Total input data size: %i" % seq_count)
    print ("Test set size: %i" % len(test_odict))
    print ("Train set size: %i" % len(train_odict))

    with open(out_fasta_test, 'w') as f:
        for keys, vals in test_odict.items():
            f.write(">" + keys + "\n" + str(vals) + "\n")
    with open(out_fasta_train, 'w') as f:
        for keys, vals in train_odict.items():
            f.write(">" + keys + "\n" + str(vals) + "\n")


def path_exists(file):
    if (not os.path.exists(file)):
        sys.exit("File path error: " + file + " does not exist!")


def args_setup():
    '''Command line arguments parsing function'''
    parser = argparse.ArgumentParser(description="Tool for splitting a fasta files into 2 separates files - train set & test set, based on a test size proportion")
    parser.add_argument("fasta", help="fasta file with nucleotide/aa sequences")
    parser.add_argument("-s", "--size_testset", metavar="<test set size (float)>", type=float, required=True,
                        help="required size (proportion, e.g. 0.2) of the test set")
    parser.add_argument("-otr", "--output_train", metavar='<path-to-output>', required=True,
                        help="path to output file name for the fasta train set")
    parser.add_argument("-ote", "--output_test", metavar='<path-to-output>', required=True,
                        help="path to output file name for the fasta test set")

    return parser.parse_args()


if __name__ == "__main__":
    main()
