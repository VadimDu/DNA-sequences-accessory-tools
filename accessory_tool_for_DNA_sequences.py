
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 19:00:39 2019
Modified on Sep 3, 2020

@author: dnx
"""
#Accessory tool for DNA sequences filtering and basic statistics; version 0.4
#Author: Dani (Vadim) Dubinsky (dani.dubinsky@gmail.com)

import argparse
import os
import sys
import csv
import collections
import numpy as np
from pathlib import Path

try:
    from Bio import SeqIO
    from Bio.SeqUtils import GC
    from Bio.SeqRecord import SeqRecord
except:
    sys.exit("This program requires Python3 Bio module, please install it")

#   
#--------------------------Main function--------------------------#
#
def main(): 
    
    args = args_setup()

    path_exists(args.fasta)
        
    seq_odict, num_seq = load_fasta(args.fasta)
    
    if args.id_file:
        path_exists(args.id_file)
        list_ids = load_id(args.id_file)
    
    if args.loc_on_contig:
        path_exists(args.loc_on_contig)
        contigs_id, coord_l, coord_h = load_contig_coords_id(args.loc_on_contig)
    
    if args.retrieve:
        try:
            assert list_ids
        except NameError:
            print("Error: IDs list file (--id-file) is not defined")
        else:
            retrieve_seq(seq_odict, list_ids, args.fasta)
      
    if args.delete:
        try:
            assert list_ids
        except NameError:
            print("Error: IDs list file (--id-file) is not defined")
        else:
            delete_seq(seq_odict, list_ids, num_seq, args.fasta)
    
    if args.trim:
        try:
            assert contigs_id
        except NameError:
            print("Error: genomic coordinates list (--loc-on-contig) is not defined")
        else:
            trim_contig(seq_odict, contigs_id, coord_l, coord_h, args.fasta)
    
    if args.min_length:
        filter_len(seq_odict, args.min_length, args.fasta)
    
    if args.basic_stats:
        basic_stats(seq_odict, args.fasta)

#
#--------------------------Utility functions--------------------------#
#

def args_setup():
    '''Command line arguments parsing function'''
    parser = argparse.ArgumentParser(description="Accessory tool for DNA sequences filtering and basic statistics, version 0.3\n By Vadim (Dani) Dubinsky (dani.dubinsky@gmail.com)")
    parser.add_argument("fasta", help = "fasta file/s with genes or contigs (supports the usage of wildcard '*' to select multiple files)", nargs = "+") #This is the positional argument (mandatory). The parser expect >=1 of positional arguments
    parser.add_argument("-i","--id-file", metavar = "<IDs list>", required = False, help = "text file with IDs list (sequence headers) - one per line") #optional arguments
    parser.add_argument("-r", "--retrieve", action = "store_true", required = False, help = "search and retrieve sequences based on IDs list") #optional arguments
    parser.add_argument("-d", "--delete", action = "store_true", required = False, help = "search and delete sequences based on IDs list") #optional arguments
    parser.add_argument("-l", "--loc-on-contig", metavar = "<genomic coordinates lists>", required = False, help = "list with target genes/contigs IDs (1st column) and specific genomic coordinates where to trim each sequence (start: 2nd column, end: 3rd column) - tab-delimeted file") #optional arguments
    parser.add_argument("-t", "--trim", action = "store_true", required = False, help = "trim specific genes/contigs by nucleotide coordinates based on an input genomic coordinates lists") #optional arguments
    parser.add_argument("-m","--min-length", metavar = "<length (bp)>" , type = int, required = False, help = "delete sequences shorter than --min-length")
    parser.add_argument("-s", "--basic-stats", action = "store_true", required = False, help = "output basic sequence statistics report (total sequences length, number of sequences, GC content, assembly N50)")
    return (parser.parse_args())

def path_exists(file):
    if type(file) == str:
        if (not os.path.exists(file)):
            #raise argparse.ArgumentTypeError("{0} does not exist".format(file))
            sys.exit("Error: " + file + " file does not exist!")
    else:
        for f in file: #file type is list, because multiple files were provided with a wildcard(*)
            if (not os.path.exists(f)):
                #raise argparse.ArgumentTypeError("{0} does not exist".format(file))
                sys.exit("Error: " + f + " file does not exist!")

def load_fasta(fasta):
    '''Load & copy the fasta file/s into an ordered dict/s
       Support several files in parallel'''
    list_of_odicts = []
    for file in fasta:
        seq_odict = collections.OrderedDict() #for storing sequences IDs (keys) and their nucleotide/amino-acid sequences (vals)
        for sequence in SeqIO.parse(file, "fasta"):
            seq_odict[str(sequence.id)] = sequence.seq
        list_of_odicts.append(seq_odict)
        num_seq = len(seq_odict)
        print("Fasta file: %s is loaded" % (os.path.split(file)[1]))
    print("Current working dir: %s" % (os.path.dirname(file))) # or os.path.split(file)[0]
    return (list_of_odicts, num_seq)

def load_id(id_file):
    '''Load a txt file with IDs list (sequence headers)'''
    with open(id_file, 'r') as f:
        list_ids = [line.strip() for line in f]
    print("IDs list file: %s is loaded" % (id_file))
    return (list_ids)

def load_contig_coords_id(loc_on_contig):
    '''Load a txt file with a table of contigs ID (1st column) and genomic coordinates (start->end, columns 2,3)'''
    contigs_id, coord_l, coord_h = [],[],[]
    with open(loc_on_contig, 'r') as f:
        reader = csv.reader(f, delimiter = "\t")
        for r in reader:
            contigs_id.append(r[0])
            coord_l.append(r[1])
            coord_h.append(r[2])
    return(contigs_id, coord_l, coord_h)

def output_fasta(retrieve, delete, min_length, seq, fasta):
    '''Copy the results to files in fasta format
       Support several files in parallel in in "min_length_mode, but only single files in "retrieve/delete" modes'''
    if (retrieve):
        with open(fasta + '_retrieved.fasta', 'w') as f:
            for keys,vals in seq.items():
                f.write(">" + keys + "\n" + str(vals) + "\n")
        print("Output file <" + os.path.split(fasta)[1] + '_retrieved.fasta' + "> was created")
    elif (delete):     
        with open(fasta + '_filtered.fasta', 'w') as f:
            for keys,vals in seq.items():
                f.write(">" + keys + "\n" + str(vals) + "\n")
        print("Output file <" + os.path.split(fasta)[1] + '_filtered.fasta' + "> was created")
    elif (min_length):
        for odict, file in zip(seq, fasta):
            with open(file + '_min_length.fasta', 'w') as f:
                for keys,vals in odict.items():
                    f.write(">" + keys + "\n" + str(vals) + "\n")
            print("Output file <" + os.path.split(file)[1] + '_min_length.fasta' + "> was created")

def retrieve_seq(seq_odict, list_ids, fasta):
    '''Search and match sequences from list_ids file (ID headers list). This method was changed in order to retrieve the sequences in the same order as in id_list file
       Works only on a single file'''
    #found_seq = {}
    found_seq = collections.OrderedDict()
    for l in list_ids:
        for keys,vals in seq_odict[0].items(): #keys are IDs, vals are the sequence itself
            #if keys in list_ids:
            if keys == l:
                found_seq[keys] = vals
    print("Search has finished")
    if not found_seq: 
        print ("No matches found!")
    else:    
        print ("Matches found! Writing into file")
        output_fasta(1,0,0, found_seq, fasta[0])

def delete_seq(seq_odict, list_ids, num_seq, fasta):
    '''Keep only the non-matched sequences (exclude the sequences from ID header list)
       Works only on a single file'''
    for key in [key for key in seq_odict[0] if key in list_ids]: del seq_odict[0][key] #Using dictionary comprehension 
    print("Search has finished")
    if len(seq_odict[0]) == num_seq:
        print ("No sequences deleted!")
    else:    
        print ("Sequences deleted! Writing into file")
        output_fasta(0,1,0, seq_odict[0], fasta[0])

def trim_contig(seq_odict, contigs_id, coord_l, coord_h, fasta):
    '''Trim specific contigs by nucleotide coordinates based on an input list with genomic coordinates (start->end)
       Works only on a single file'''
    seqs = []
    for keys,vals in seq_odict[0].items():
        if (keys not in contigs_id): #if this contig is not in the input list of contigs to trim, then copy it without changes
            seqs.append(SeqRecord(vals, id=str(keys), description=''))
        for i in range(len(contigs_id)):
            if (keys == contigs_id[i]):
                seqs.append(SeqRecord(vals[0:int(coord_l[i])-1] + vals[int(coord_h[i])::], id=str(keys), description=''))
                     
    SeqIO.write(seqs, fasta[0] + "_trimmed.fasta", "fasta")
    print("Output file <" + os.path.split(fasta[0])[1] + '_trimmed.fasta' + "> was created")

def filter_len(seq_odict, min_len, fasta):
    '''Delete sequences shorter than min_len
       Works on several files in parallel (with * wildcard)'''
    for odict in seq_odict:
        seq_del = []
        for keys,vals in odict.items():
            if len(vals) < (min_len):
                seq_del.append(keys)
        for i in seq_del:
            del odict[i]
    output_fasta(0,0,1, seq_odict, fasta)

def basic_stats(seq_odict, fasta):
    '''Calculation of basic sequence statistics, including N50 of contigs. The N50 is defined as the minimum contig length needed to cover 50% of the genome
       Works on several files in parallel (with * wildcard)'''
    list_length, list_num, list_gc, list_n50, list_id = [],[],[],[],[]
    for odict, file in zip(seq_odict, fasta):
        seq_len, seq_gc = [],[]
        for i in odict: #loop over all the sequences in a fasta file (=odict)
            seq_len.append(len(odict[i]))
            seq_gc.append(GC(odict[i]))
        list_id.append(Path(file).stem)
        list_length.append(round(sum(seq_len) / 1e6,3))
        list_num.append(len(seq_len))
        list_gc.append(round(np.median(seq_gc),2))
        
        #N50 calculation
        seq_len_sorted = sorted(seq_len, reverse=True) # sort contigs longest>shortest
        csum = np.cumsum(seq_len_sorted) #Cumulative sum of the lengths of all the contigs
        n_div_2 = int(sum(seq_len_sorted)/2) #The 50% of the total assembly length 
        csum_n_div_2 = min(csum[csum >= n_div_2]) #The cumulative sum of the lengths of the assembly that contain 50% of the total assembly length
        ind = np.where(csum == csum_n_div_2) #The index for cumsum >= N/2 (50% of the total assembly length)
        list_n50.append(seq_len_sorted[int(ind[0])]/1000) #to convert to Kb

    #Using csv.writer to output the results to a csv file:
    os.chdir(os.path.split(fasta[0])[0])
    with open("Bio-stats_report.txt", "w") as f:
        colnames = ["Fasta file with genes/contigs","Gene/contig length [Mb]","Number of genes/contigs","Median GC content", "N50 [Mb]"]
        writer = csv.writer(f, delimiter = "\t") #sys.stderr 
        writer.writerow(colnames)
        for i in range(len(list_length)):
            writer.writerow([list_id[i], list_length[i], list_num[i], list_gc[i], list_n50[i]])
    print("<Bio-stats_report.txt> with sequence bio-statistics was created in: " + os.path.split(fasta[0])[0])

#--------------------------END--------------------------------------#      

if __name__ == "__main__":
    main()
