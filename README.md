# DNA sequences accessory tool
An accessory for performing basic manipulations and obtaining bio-statistical information on DNA sequences. Can work also on amino-acid sequences.

## Python modules requirements
You need to have Python version >=3.0 and the following modules installed:
<br/>numpy
<br/>Bio

## Usage instructions
Run the provided 'accessory_tool_for_DNA_sequences.py' script with Python3 on a fasta file (or on multiple fasta files using a wildcard) that contain either nucleotide or amino-acid sequences in standard fasta format. 
With this script you can make summary bio-statistics (number of sequences, total sequences length, GC content and N50 value). You can retrieve or delete specific sequences based on sequence ID (headers) list. You can filter the fasta file to delete sequences shorter than a defined threshold (X bp). Finally, you can trim specific sequences based genomic coordinates (start position, end position).
