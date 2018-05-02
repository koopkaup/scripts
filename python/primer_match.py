#!/usr/bin/python3

from __future__ import nested_scopes

'''
Run as:
primer_match.py sisendfail (fasta järjestused) fprimer rprimer mismatch error both.out forw.out rev.out primer.out
'''
__author__ = "Kristjan Oopkaup, Jens-Konrad Preem"
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Kristjan Oopkaup"
__email__ = "kristjan.oopkaup@gmail.com"
__status__ = "Production"

import sys
import string
import time
from time import localtime, strftime

import regex
try:
    from Bio import SeqIO
except ImportError:
    raise ImportError("BioPython is not installed!")


#####  Functions  #####

def multiple_replace(dict, text): 
    '''Teeb stringis asendused vastavalt antud asendus-raamatukogule'''
    # Create a regular expression  from the dictionary keys
    dictRegex = regex.compile("(%s)" % "|".join(map(regex.escape, dict.keys())))
    # For each match, look-up corresponding value in dictionary
    return dictRegex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text) 

def rc(dna):
    '''p88rdkomplemendi meetod'''
    complements = str.maketrans('acgtrymkbdhvACGTRYMKBDHV', 'tgcayrkmvhdbTGCAYRKMVHDB')
    rcseq = dna.translate(complements)[::-1]
    return rcseq

#####  Program  #####

start = time.time()
print('Primer search started at: %s' % (strftime("%H:%M:%S", localtime())))

infile=sys.argv[1] #fasta järjestused
fprim=sys.argv[2] #forward primer
rprim=sys.argv[3] #reverse primer
mismatch=sys.argv[4] #mismatch number
errors=sys.argv[5] #error type: insertion, deletion, substitution, everything (i, d, s, e)
outfile1=sys.argv[6] #both primers match
outfile2=sys.argv[7] #forward primers match
outfile3=sys.argv[8] #reverse primers match
outfile4=sys.argv[9] #primers locations

f = open(infile, 'r')
dout = open(outfile1, 'w+')
fout = open(outfile2, 'w+')
rout = open(outfile3, 'w+')
ploc = open(outfile4, 'w+')

#seqlist = list(SeqIO.parse(f, "fasta"))
#pikkus=len(seqlist)

#rprim p88rdkomplement

rprim=rc(rprim)

#k6du asenduste raamatukogu
dict = {
"R" : "[GA]",
"Y" : "[TC]",
"M" : "[AC]",
"K" : "[GT]",
"S" : "[GC]",
"W" : "[AT]",
"H" : "[ACT]",
"B" : "[GTC]",
"V" : "[GCA]",
"D" : "[GTA]",
"N" : "[GTAC]",
} 

#asendused
fprima=multiple_replace(dict, fprim.upper()) 
rprima=multiple_replace(dict, rprim.upper())
#fuzzy regex praimeritest j2rjestustele sobitamiseks
fex = regex.compile('('+fprima+'){'+errors+'<='+mismatch+'}')
rex = regex.compile('('+rprima+'){'+errors+'<='+mismatch+'}')
#fex = regex.compile(r'(^'+fprima+'|'+fprima+'$){'+errors+'<='+mismatch+'}')
#rex = regex.compile(r'(^'+rprima+'|'+rprima+'$){'+errors+'<='+mismatch+'}')

#sobitamine ja vastavasse v2ljundfaili print
bcount = 0
fcount = 0
rcount = 0

ploc.write("**** Primers summary ****\n")
ploc.write("Forward primer: %s (%i nt)\n" % (fprim, len(fprim)))
ploc.write("Reverse primer: %s (%i nt)\n\n" % (rprim, len(rprim)))
for record in SeqIO.parse(f, "fasta"):
    #m=regex.findall(fex, str(record.seq), overlapped=True)
    #n=regex.findall(rex, str(record.seq), overlapped=True)
    m = regex.search(fex, str(record.seq), overlapped=True)
    n = regex.search(rex, str(record.seq), overlapped=True)
    if m and n:
        SeqIO.write(record, dout, "fasta")
        bcount += 1
        ploc.write("*"*24 + '\n')
        ploc.write( "%s (length: %i bp)\n" % (record.id, len(record.seq)))
        ploc.write("Forward primer match:\n")
        ploc.write("\n")
        for matchf in regex.finditer(fex, str(record.seq), overlapped=True):
            ploc.write("Found: %s (%i nt)\n" % (matchf.group(0),len(matchf.group(0))))
            ploc.write("At location: from %i nt to %i nt\n" % (matchf.start()+1, matchf.end()+1))
            ploc.write("\n")
            ploc.write("-"*24 + '\n')
            ploc.write("Reverse primer match:\n")
            ploc.write("\n")
            for matchr in regex.finditer(rex, str(record.seq), overlapped=True):
                ploc.write("Found: %s (%i nt)\n" % (matchr.group(0),len(matchr.group(0))))
                ploc.write("At location: from %i nt to %i nt\n" % (matchr.start()+1, matchr.end()+1))
                ploc.write("\n")
                if m:
                    SeqIO.write(record, fout, "fasta")
                    fcount += 1
                    if n:
                        SeqIO.write(record, rout, "fasta")
                        rcount += 1

                        print("Filename: %s\n" % (infile.split('/')[-1] if '/' in infile else infile))
#print("Total number of reads: " + str(pikkus) + "\n")
print("Mismatches allowed: " + str(mismatch) + "\n")
print("Both primers match: " + str(bcount) + "\n")
print("Forward primers match: " + str(fcount) + "\n")
print("Reverse primers match: " + str(rcount) + "\n")

dout.close()
fout.close()
rout.close()
ploc.close()
f.close()

end = time.time()
print("Time elapsed: %i minutes" % (round((end - start)/60)))