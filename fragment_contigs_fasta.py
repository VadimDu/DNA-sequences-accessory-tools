#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  2 13:11:57 2021

@author: dnx
"""
# title: fragment_contigs_fasta.py
# description: Tool for fragmentation of fasta files (nucleotide contigs in this case) into equal length pieces
# author: Vadim (Dani) Dubinsky [dani.dubinsky@gmail.com]
# date: 02/05/2021
# required modules: Biopython

import argparse
import os
import sys
from collections import OrderedDict

try:
    from Bio import SeqIO
except:
    sys.exit("This program requires Python3 Bio module, please install it")


def main():
    
    args = args_setup()
    
    path_exists(args.fasta)
    
    #Constant
    FRAG_LEN = 10000
    
    if args.fragment_len:
        FRAG_LEN = args.fragment_len
    
    fragment_fasta(args.fasta, args.output, FRAG_LEN)


def fragment_fasta(in_fasta, out_fasta, frag_len):
    ''' Fragment contigs (sequences) in a fasta file into equal length fragments of a defined length
    Parameters
    ----------
    in_fasta : input fasta file (path)
    out_fasta : output file for the fragmented contigs (path)
    frag_len : required fragment length to split the input fasta file to (this equally sized fragments)
    Returns
    -------
    None.
    '''
    input_contigs_len = 0
    input_contigs_count = 0
    frag_total_len = 0
    contig_frags = OrderedDict()
    
    for contig in SeqIO.parse(in_fasta, "fasta"):
        input_contigs_len += len(contig.seq)
        input_contigs_count += 1
        if len(contig.seq) < frag_len: continue #skip contigs shorter than the selected fragment length
    
        frag_count = 0
        for idx in range(0, len(contig.seq), frag_len):
            fragment = contig.seq[idx:(idx+frag_len)]
            if len(fragment) < frag_len: continue ##skip contig fragments shorter than the selected fragment length
            frag_count += 1
            frag_total_len += len(fragment)
            contig_frags[str(contig.description) + "_" + str(frag_count)] = fragment #keep the original sequence header and append contig counter
    
    print("Original input contigs number and length: {}, {} bp \nFragmented contigs number and length: {}, {} bp".format(input_contigs_count, input_contigs_len, len(contig_frags), frag_total_len))
    
    with open(out_fasta, 'w') as f:
        for keys,vals in contig_frags.items():
            f.write(">" + keys + "\n" + str(vals) + "\n")


def path_exists(file):
        if (not os.path.exists(file)):
            sys.exit("File path error: " + file + " does not exist!")


def args_setup():
    '''Command line arguments parsing function'''
    parser = argparse.ArgumentParser(description="Tool for fragmentation of fasta files (nucleotide contigs in this case) into equal length pieces (By Vadim-Dani Dubinsky dani.dubinsky@gmail.com)")
    parser.add_argument("fasta", help = "fasta file with nucleotide contigs (also accepts amino-acid sequences)")
    parser.add_argument("-f","--fragment_len", metavar = "<fragment length (bp)>" , type = int, required = False, default = 10000, help = "the required fragment length to split the input fasta file to equally sized fragments (10kb by default)")
    parser.add_argument("-o", "--output", metavar = '<path-to-output>', required = True, help = "path to output file name for the fragmented contigs")
    return (parser.parse_args())


if __name__ == "__main__":
    main()