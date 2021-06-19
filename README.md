# DNA sequences accessory tools
An accessory tool for performing basic manipulations and obtaining bio-statistical information on DNA sequences. Can work also on amino-acid sequences:
`accessory_tool_for_DNA_sequences.py`


A tool for fragmentation of fasta files (usefull in case of nucleotide contigs) into equal length fragments:
`fragment_contigs_fasta.py`

## Python modules requirements
You need to have Python version >=3.0 and the following modules installed:
<br/>numpy
<br/>Bio

## Usage instructions
Run the provided `accessory_tool_for_DNA_sequences.py` script with Python3 on a fasta file (or on multiple fasta files using a wildcard) that contain either nucleotide or amino-acid sequences in standard fasta format.<br/> 
With this script you can make summary bio-statistics (number of sequences, total sequences length, GC content and N50 value). You can retrieve or delete specific sequences based on sequence ID (headers) list. You can filter the fasta file to delete sequences shorter than a defined threshold (X bp). Finally, you can trim specific sequences (remove nucleotides) based on genomic coordinates (start position, end position over each sequence).<br/>

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
                        list with target genes/contigs IDs (1st column) and specific genomic coordinates where to trim each sequence (start: 2nd column, end: 3rd column) - tab-delimeted file
  -t, --trim            trim specific genes/contigs by nucleotide coordinates based on an input genomic coordinates lists
  -m <length (bp)>, --min-length <length (bp)>
                        delete sequences shorter than --min-length
  -s, --basic-stats     output basic sequence statistics report (total sequences length, number of sequences, GC content, assembly N50)
```
**Note:** the following options can work on multiple fasta files provided with a wildcard: --min-length and --basic-stats. The rest of the options work only on a single fasta file.<br/>

## Examples
* To delete sequence shorter than 1000bp from a fasta file:<br/>
`python3 accessory_tool_for_DNA_sequences.py your_file.fasta -m 1000`<br/>
Output file <your_file.fasta_min_length.fasta> was created<br/>
* To retrieve 3 specific sequences from your fasta file, provide a text file with these 3 seqeuence headers:<br/>
```
cat header_ids.txt
NODE_1_length_34583
NODE_3_length_33269
NODE_5_length_28153

python3 accessory_tool_for_DNA_sequences.py your_file.fasta --id-file header_ids.txt --retrieve
Output file <your_file.fasta_retrieved.fasta> was created
```
* To trim 4 specific sequences according to genomic coordinates (remove nucleotides at specific positions in each sequence, start pos. -> end pos.). Provide a tab-delimeted file with 3 columns: sequence ID, start position, end position:<br/>

```
cat sequence_id_coords.txt
NODE_1_length_53121	10	20
NODE_2_length_52466	15	25
NODE_3_length_50324	25	40
NODE_10_length_43224	1	50

python3 accessory_tool_for_DNA_sequences.py your_file.fasta --loc-on-contig sequence_id_coords.txt --trim
Output file <your_file.fasta_trimmed.fasta> was created
```

## Fragment a fasta file (with nucleotide contigs) into equal length fragments (user input)
```
usage: fragment_contigs_fasta.py [-h] [-f <fragment length bp>] -o <path-to-output> fasta

Tool for fragmentation of fasta files (nucleotide contigs in this case) into equal length pieces (By Vadim-Dani Dubinsky dani.dubinsky@gmail.com)

positional arguments:
  fasta                 fasta file with nucleotide contigs (also accepts amino-acid sequences)

optional arguments:
  -h, --help            show this help message and exit
  -f <fragment length (bp)>, --fragment_len <fragment length (bp)>
                        the required fragment length to split the input fasta file to equally sized fragments (10kb by default)
  -o <path-to-output>, --output <path-to-output>
                        path to output file name for the fragmented contigs
```

## Random sampling (without replacment) of fasta files, required number of sequences by user input
```
usage: random_sampling_fasta.py [-h] -s <sequence number int> -o <path-to-output> fasta

Tool for random sampling (without replacment) of fasta files

positional arguments:
  fasta                 fasta file with nucleotide/aa sequences

optional arguments:
  -h, --help            show this help message and exit
  -s <sequence number (int)>, --sequence_num <sequence number (int)>
                        the required number of sequences to sample from the input fasta
  -o <path-to-output>, --output <path-to-output>
                        path to output file name for the sampled subset of the original fasta
```

## Split a fasta files into 2 separates files - train set & test set, based on a test size proportion (user input)
```
usage: split_train_test_fasta.py [-h] -s <test set size float> -otr <path-to-output> -ote <path-to-output> fasta

Tool for splitting a fasta files into 2 separates files - train set & test set, based on a test size proportion

positional arguments:
  fasta                 fasta file with nucleotide/aa sequences

optional arguments:
  -h, --help            show this help message and exit
  -s <test set size (float)>, --size_testset <test set size (float)>
                        required size (proportion, e.g. 0.2) of the test set
  -otr <path-to-output>, --output_train <path-to-output>
                        path to output file name for the fasta train set
  -ote <path-to-output>, --output_test <path-to-output>
                        path to output file name for the fasta test set
```

## Handling of none-ATGC nucleotides in a fasta file by replacing them with a random nucleotides
```
usage: handle_none_ATGC.py [-h] -o <path-to-output> fasta

Tool for handling none-ATGC nucleotides in a fasta file by replacing them with a random nucleotides

positional arguments:
  fasta                 fasta file with nucleotide contigs (also accepts amino-acid sequences)

optional arguments:
  -h, --help            show this help message and exit
  -o <path-to-output>, --output <path-to-output>
                        path to output file name for the corrected non-ATGC fasta
```
