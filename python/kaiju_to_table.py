#!/usr/bin/env python3

'''
Merge kaiju summary report files for a given taxonomic rank to single table. Can be used to separate bacteria and archaea.
Proportions are calculated separately for bacteria and archaea.
kaijuReport must be run using argument -p (Print full taxon path)
for example: kaijuReport -t nodes.dmp -n names.dmp -r phylum -p -i kaiju.out

Run as:
kaiju_to_table.py input_files
'''
__author__ = "Kristjan Oopkaup"
__license__ = "GPL"
__version__ = "0.7"
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
    elif kingdom == 'Both':
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

kingdom = int(input('Do you want bacteria [1], archaea [2] or both [3]?\n'))
taxonLvl = int(input('What taxonomic rank you have in your input files?\n'
    'Phylum [1],Class [2],Order [3],Family [4],Genus [5],Species [6]\n'))
numbers = int(input('Do you want to represent proportions [1] or sequence counts [2]?\n'))

taxa = ['kingdom','phylum','class','order','family','genus','species']

if kingdom == 1 and numbers == 1:
    table_out = open('kaiju_table_bacteria_proportions_'+taxa[taxonLvl]+'.csv', 'w+')
elif kingdom == 1 and numbers == 2:
    table_out = open('kaiju_table_bacteria_count_'+taxa[taxonLvl]+'.csv', 'w+')
elif kingdom == 2 and numbers == 1:
    table_out = open('kaiju_table_archaea_proportions_'+taxa[taxonLvl]+'.csv', 'w+')
elif kingdom == 2 and numbers == 2:
    table_out = open('kaiju_table_archaea_count_'+taxa[taxonLvl]+'.csv', 'w+')
elif kingdom == 3 and numbers == 1:
    table_out = open('kaiju_table_proportions_'+taxa[taxonLvl]+'.csv', 'w+')
elif kingdom == 3 and numbers == 2:
    table_out = open('kaiju_table_count_'+taxa[taxonLvl]+'.csv', 'w+')


taxaDict = defaultdict(list)
firstLine = taxa[taxonLvl]
totalSum = list()
for i in range(1, len(sys.argv)):
    with open(sys.argv[i]) as kaijuIN:
        total = 0
        firstLine += '\t%s' % (sys.argv[i])
        lines = kaijuIN.readlines()[2:-5]
        for line in lines:
            if kingdom == 1:
                taxaDict, total = taxaFunc(line, 'Bacteria', total, i)
            elif kingdom == 2:
                taxaDict, total = taxaFunc(line, 'Archaea', total, i)
            elif kingdom == 3:
                taxaDict, total = taxaFunc(line, 'Both', total, i)
        for k, v in taxaDict.items():
            for zeros in range(1, i-len(v)+1):
                taxaDict[k].append(0)
        totalSum.append(total)
        if numbers == 1:
            for k, v in taxaDict.items():
                v = v[:-1] + [round(v[-1]/total,5)]
                taxaDict[k] = v

table_out.write(firstLine + '\n')
for taxon, nr in taxaDict.items():
    table_out.write('%s\t%s\n' % (taxon, '\t'.join(str(v) for v in nr)))
table_out.write('Sequence count\t%s' % ('\t'.join(str(v) for v in totalSum)))

table_out.close()