#!/usr/bin/python3

'''
Extract blast search hits from query fasta file and write a new fasta file.
Run as:
blast_to_fasta.py blast.results.csv query.fasta output.fasta
'''
__author__ = "Kristjan Oopkaup"
__license__ = "GPL"
__version__ = "0.1"
__email__ = "kristjan.oopkaup@gmail.com"
__status__ = "Production"

import sys
try:
    from Bio import SeqIO, SearchIO
except ImportError:
    raise ImportError("BioPython is not installed!")

blastfasta = open(sys.argv[3], 'w+') # Output fasta file

query_id = set()

for i in SearchIO.parse(sys.argv[1], 'blast-tab'): # Blast results file
	query_id.add(i.id)

for y in SeqIO.parse(sys.argv[2], 'fasta'): # Query fasta
	if y.id in query_id:
		SeqIO.write(y, blastfasta, 'fasta')

blastfasta.close()