#!/usr/bin/python
"""
Merges expression values for overlapping group of genes. See for details:

Christelle, and Watson. "Errors in RNA-Seq quantification affect genes of
relevance to human disease." Genome Biology 16.1 (2015): 177.

# List mmquant.in:

gene	sample1	sample2
gene00	1	2
gene02	1	2
gene05	1	2
gene07	10	20
gene08	10	20
gene09	10	20
gene04	10	20
gene04_gene41_gene06_gene48_gene43_gene31_gene68	1	2
gene04_gene41_gene48_gene43_gene68	1	2
gene04_gene41_gene48_gene43_gene68_gene84	1	2
gene04_gene48_gene77_gene43_gene68_gene58	1	2
gene04_gene48_gene53_gene43_gene68	1	2
gene04_gene48_gene43_gene68	1	2
gene04_gene48_gene43_gene68_gene44	10	20
gene04_gene84	0	0
gene05	0	0
gene41	10	20
gene53	10	20
gene53_gene88	10	20
gene53_gene39	10	20
gene53_gene39_gene88	1	20
gene53_gene39_gene88_gene84	1	0
gene53_gene39_gene84	10	20
gene53_gene86	1	20
gene53_gene86_gene39	1	0
gene53_gene86_gene39_gene88	10	20
gene53_gene86_gene39_gene88_gene84	10	20
gene53_gene86_gene39_gene84	10	20


# output mmquant.out:

gene	sample1	sample2
gene00	1	2
gene02	1	2
gene05	1	2
gene07	10	20
gene08	10	20
gene09	10	20
gene04	10	20
gene05	0	0
gene41	10	20
gene53	10	20
gene04_gene48_gene43_gene68_gene44	10	20
gene53_gene86_gene39_gene88_gene84	62	160

# command:

python mergeGeneGroups_mmquant.py -i mmquant.in -o mmquant.out -l 2 -c 10

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import argparse, sys  # for input options
import collections

############################# options #############################

class CommandLineParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)  # command line syntax errors


parser = CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the input file',
    type=str,
    required=True)
parser.add_argument(
    '-o',
    '--output',
    help='name of the output file',
    type=str,
    required=True)
parser.add_argument(
    '-c',
    '--count',
    help='minimum number of reads per group to keep it (optional)',
    default=0,
    type=int,
    required=False)
parser.add_argument(
    '-l',
    '--length',
    help=
    'minimum number of genes per group to merge it with larger group (optional)',
    default=1,
    type=int,
    required=False)
args = parser.parse_args()

############################# functions ###########################

def mergeMMG(Dict):
    DictKeys = sorted(Dict.keys(), key=len, reverse=True)
    dictNew = Dict.copy()
    for m in DictKeys:
        # print()
        overlap = 0
        for n in DictKeys:
            mm = set(m.split('_'))
            nn = set(n.split('_'))
            # print('m n:', m, n)
            if mm.issubset(nn) and m != n and n in dictNew:
                # print('m n:', m, n)
                key2rm = m
                key2keep = n
                overlap += 1
        if overlap == 1:
            newVal = []
            for i, j in zip(Dict[key2rm], Dict[key2keep]):
                newVal.append(sum([i, j]))
            # print(key2rm, Dict[key2rm], key2keep, Dict[key2keep])
            Dict[key2keep] = newVal
            dictNew[key2keep] = newVal
            del dictNew[key2rm]
    return dictNew

def list2print(Lits):
    LitsP = '\t'.join(str(e) for e in Lits)
    return LitsP

############################# program #############################

output = open(args.output, 'w')
exprDic = {}

with open(args.input) as datafile:
    header = datafile.readline()
    output.write(header)

    for line in datafile:
        words = line.split()
        genes = words[0]
        counts = list(map(int, words[1:]))
        genesNumber = len(genes.split('_'))
        if ('_' in genes and sum(counts) >= args.count) or '_' not in genes:
            if genesNumber >= args.length:
                exprDic[genes] = counts
            else:
                countsP = list2print(counts)
                output.write("%s\t%s\n" % (genes, countsP))

    exprDicMerged = mergeMMG(exprDic)
    for key, values in exprDicMerged.items():
        valuesP = list2print(values)
        output.write("%s\t%s\n" % (key, valuesP))

datafile.close()
output.close()
print('Done!')
