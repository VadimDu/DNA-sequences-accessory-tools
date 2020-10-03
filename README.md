# DNA sequences accessory tool
An accessory for performing basic manipulations and obtaining bio-statistical information on DNA sequences. Can work also on amino-acid sequences.

## Python modules requirements
You need to have Python version >=3.0 and the following modules installed:
<br/>numpy
<br/>Bio

## Usage instructions
Run the provided `accessory_tool_for_DNA_sequences.py` script with Python3 on a fasta file (or on multiple fasta files using a wildcard) that contain either nucleotide or amino-acid sequences in standard fasta format.<br/> 
With this script you can make summary bio-statistics (number of sequences, total sequences length, GC content and N50 value). You can retrieve or delete specific sequences based on sequence ID (headers) list. You can filter the fasta file to delete sequences shorter than a defined threshold (X bp). Finally, you can trim specific sequences based genomic coordinates (start position, end position).

## Full command-line options (--help)
```
usage: accessory_tool_for_DNA_sequences.py [-h] [-i <IDs list>] [-r] [-d] [-l <genomic coordinates lists>] [-t] [-m <length bp>] [-s] fasta [fasta ...]

Accessory tool for DNA sequences filtering and basic statistics, version 0.3 By Vadim (Dani) Dubinsky (dani.dubinsky@gmail.com)

positional arguments:
  fasta                 fasta file/s with genes or contigs (supports the usage of wildcard '*' to select multiple files)

optional arguments:
  -h, --help            show this help message and exit
  -i <IDs list>, --id-file <IDs list>
                        text file with IDs list (sequence headers) - one per line
  -r, --retrieve        search and retrieve sequences based on IDs list
  -d, --delete          search and delete sequences based on IDs list
  -l <genomic coordinates lists>, --loc-on-contig <genomic coordinates lists>
                        list with target contigs IDs (1st column) and specific genomic coordinates where to trim each contig (start: 2nd column, end: 3rd column) - tab-delimeted file
  -t, --trim            trim specific contigs by nucleotide coordinates based on an input genomic coordinates lists
  -m <length (bp)>, --min-length <length (bp)>
                        delete sequences shorter than --min-length
  -s, --basic-stats     print basic sequence statistics (total sequences length, number of sequences, GC content, assembly N50)
```

## Examples
To delete sequence shorter than 1000bp from a fasta file:<br/>
`Python3 accessory_tool_for_DNA_sequences.py your_file.fasta -m 1000`
Output file <ERS2664665_scaffolds_mod.fa_min_length.fasta> was created
<br/>
To retrieve 3 specific sequences from your fasta file, provide a text file with these 3 seqeuence headers:
```
cat header_ids.txt
NODE_1_length_34583
NODE_3_length_33269
NODE_5_length_28153

Python3 accessory_tool_for_DNA_sequences.py your_file.fasta --id-file header_ids.txt --retrieve
Output file <your_file.fasta_retrieved.fasta> was created
```






`
