#!/usr/bin/python

"""
Annotates the results of the differential expression analysis of genes and 
multi-mapped groups (MMG).

For MMG see:
Christelle, and Watson. "Errors in RNA-Seq quantification affect genes of
relevance to human disease." Genome Biology 16.1 (2015): 177.

# DE_mmg.csv:

genes	log2FoldChange	padj
ENSCAFG00000006589	-0.925421076069308	2.76491413073952E-06
ENSCAFG00000028960_ENSCAFG00000007572	4.92712501133988	2.77710999466866E-06
ENSCAFG00000028653_ENSCAFG00000032286_ENSCAFG00000031239_ENSCAFG00000030588_ENSCAFG00000032684	-4.45847235603138	3.08644735967747E-06


# annotation.csv:

genes	mean_pvalues	mean_pvalues_100K	mean_pbs	mean_pbs_100K
ENSCAFG00000000001	0.578429052279	0.570328779485	-0.00978219087713	-0.01044646188
ENSCAFG00000000002	0.572798116061	0.529718584861	-0.0136947483174	-0.01146510902
ENSCAFG00000007572	0.556622049299	0.661606216734	-0.0960655826936	-0.0746169332648
ENSCAFG00000028653	0.01789816	0.01398802	Inf	1.060251
ENSCAFG00000028960	0.5871494	0.6034695	-0.02696316	-0.03729555

# output.csv:

genes	log2FoldChange	padj	mean_pvalues	mean_pvalues_100K	mean_pbs	mean_pbs_100K	min_mean_pvalues	min_mean_pvalues_100K	max_mean_pbs	max_mean_pbs_100K
ENSCAFG00000006589	-0.925421076069308	2.76491413073952E-06	NA	NA	NA	NA	NA	NA	NA	NA
ENSCAFG00000028960_ENSCAFG00000007572	4.92712501133988	2.77710999466866E-06	0.5871494;0.556622049299	0.6034695;0.661606216734	-0.02696316;-0.0960655826936	-0.03729555;-0.0746169332648	0.556622049299	0.6034695	-0.02696316	-0.03729555
ENSCAFG00000028653_ENSCAFG00000032286_ENSCAFG00000031239_ENSCAFG00000030588_ENSCAFG00000032684	-4.45847235603138	3.08644735967747E-06	0.01789816;NA;NA;NA;NA	0.01398802;NA;NA;NA;NA	Inf;NA;NA;NA;NA	1.060251;NA;NA;NA;NA	0.01789816	0.01398802	NA	1.060251

# command:

$ python annotate_mmquant_DEresults.py \
    -i DE_mmg.csv \
    -a annotation.csv \
    --min mean_pvalues,mean_pvalues_100K \
    --max mean_pbs,mean_pbs_100K \
    -o output.csv

# contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

"""

############################# modules #############################

import argparse, sys  # for input options

############################# options #############################

class CommandLineParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

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
    '-a',
    '--annotation',
    help='annotation for genes',
    type=str,
    required=True)
parser.add_argument(
    '-s',
    '--min',
    help='which columns to use to extract minimum number',
    type=str,
    required=False)
parser.add_argument(
    '-l',
    '--max',
    help='which columns to use to extract maximum number',
    type=str,
    required=False)
args = parser.parse_args()


try:
    minCol = args.min.split(",")
except:
    pass

try:
    maxCol = args.max.split(",")
except:
    pass

############################# functions ###########################

def list2print(Lits):
    LitsP = '\t'.join(str(e) for e in Lits)
    return LitsP

def indexCol(colNames, header_words):
    colIndex = []
    for i in colNames:
        indnumber = header_words.index(i)
        colIndex.append(indnumber)
    return colIndex

def headerMinMax(minmax, header):
    if minmax+'Col' in globals():
        minmaxCol = eval(minmax+'Col')
        for c in minmaxCol:
            header.append(minmax+"_"+c)
    return(header)

def colMinMax(minmax, header, annot):
    if minmax+'Col' in globals():
        minmaxCol = eval(minmax+'Col')
        colIndex = indexCol(minmaxCol, header[1:])
        for s in range(len(colIndex)):
            colList = annot[colIndex[s]].split(";")
            colListNoNA = [i for i in colList if str(i) not in ['NA','Inf']]
            try:
                annot.append(eval(minmax)(list(map(float, colListNoNA))))
            except:
                annot.append('NA')
    return(annot)

############################# program #############################

output = open(args.output, 'w')

with open(args.annotation) as annotFile:
    annotHeader = annotFile.readline().rstrip().split('\t')

    annotHeader = headerMinMax('min', annotHeader)
    annotHeader = headerMinMax('max', annotHeader)
    annotHeaderP = list2print(annotHeader[1:])

    annotDic = {}
    for line in annotFile:
        annotWords = line.rstrip().split('\t')
        annotGene = annotWords[0]
        annot = annotWords[1:]
        annotDic[annotGene] = annot
annotFile.close()

with open(args.input) as deFile:
    deHeader = deFile.readline().rstrip()
    output.write("%s\t%s\n" % (deHeader, annotHeaderP))

    for line in deFile:
        words = line.rstrip().split('\t')
        mmg = words[0]
        genes = mmg.split('_')
        stats = words[1:]
        
        for g in genes:
            if 'mmgAnnot' in globals():
                for i in range(len(mmgAnnot)):
                    try:
                        addAnnot = annotDic[g][i][:]
                    except:
                        addAnnot = 'NA'
                    mmgAnnot[i] = str(mmgAnnot[i]) + ";" + str(addAnnot)
            else:
                try:
                    mmgAnnot = annotDic[g][:]
                except:
                    mmgAnnot = ['NA']*len(annot)
            # print(genes, g, mmgAnnot[0])

        mmgAnnot = colMinMax('min', annotHeader, mmgAnnot)
        mmgAnnot = colMinMax('max', annotHeader, mmgAnnot)

        mmgAnnotP = list2print(mmgAnnot)
        statsP = list2print(stats)
        output.write("%s\t%s\t%s\n" % (mmg, statsP, mmgAnnotP))
        del mmgAnnot

deFile.close()
output.close()
print('Done!')
