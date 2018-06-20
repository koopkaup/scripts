#!/usr/bin/env python3

'''
Run as:
kaiju_to_table.py input_files
'''
__author__ = "Kristjan Oopkaup"
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Kristjan Oopkaup"
__email__ = "kristjan.oopkaup@gmail.com"
__status__ = "Production"

from collections import defaultdict
import sys

####  Functions  ####

def taxaFunc(line, kingdom, total, i):
    global taxaDict
    if line.split('\t')[2].split(';')[1].strip() == kingdom:
        taxon = line.split('\t')[2].split(';')[-2].strip()
        absolute = int(line.split('\t')[1])
        if taxon not in taxaDict:
            for zeros in range(1, i):
                taxaDict[taxon].append(0)
            taxaDict[taxon].append(absolute)
        else:
            taxaDict[taxon].append(absolute)
        total += absolute
    else:
        pass
    return taxaDict, total

####  Program  ####

kingdom = int(input('Do you want bacteria [1] or archaea [2]?\n'))
taxonLvl = int(input('What taxonomy level you have in your input files?\n'
    'Phylum [1],Class [2],Order [3],Family [4],Genus [5],Species [6]\n'))
numbers = int(input('Do you want to represent proportions [1] or absolute numbers [2]?\n'))

taxa = ['Kingdom','Phylum','Class','Order','Family','Genus','Species']

if kingdom == 1:
    table_out = open('kaiju_table_bacteria_'+taxa[taxonLvl]+'.csv', 'w+')
elif kingdom == 2:
    table_out = open('kaiju_table_archaea_'+taxa[taxonLvl]+'.csv', 'w+')

taxaDict = defaultdict(list)
firstLine = taxa[taxonLvl]
for i in range(1, len(sys.argv)):
    with open(sys.argv[i]) as kaijuIN:
        total = 0
        firstLine += '\t%s' % (sys.argv[i])
        lines = kaijuIN.readlines()[2:-5]
        #lines = kaijuIN.readlines()
        #target = lines.index(-------------------------------------------)
        #print(target)
        for line in lines:
            if kingdom == 1:
                taxaDict, total = taxaFunc(line, 'Bacteria', total, i)
            elif kingdom == 2:
                taxaDict, total = taxaFunc(line, 'Archaea', total, i)
        if numbers == 1:
            for v in taxaDict.values():
                v[-1] = round(v[-1] / total,5)
            taxaDict = {k: v for k, v in taxaDict.items()}


table_out.write(firstLine + '\n')
for taxon, nr in taxaDict.items():
    if len(nr) < i:
        for zeros in range(1, i-len(nr)+1):
            taxaDict[taxon].append(0)
    table_out.write('%s\t%s\n' % (taxon, '\t'.join(str(v) for v in nr)))

table_out.close()