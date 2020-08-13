#!/usr/bin/env python

"""

script used to attribute to each genomic position the tolerance scores.

use : input | python placeHeptamerToleranceScoreFromDicotoGenomicCoordinate_onlyACGT.py dico.txt > output

where dico.txt is a dictionary file containing the heptamers (column1) and the associated tolerance score (column2)

"""


import sys
import itertools

iBED=sys.stdin
oBED=sys.stdout
scores=open(sys.argv[1], 'r')
k=int(sys.argv[2])

kmerdict={}  # creates a dictionnary that will contain all possible kmers
nct='ACTG'   # nct stands for nucleotide

for kmer in list(itertools.product(nct, repeat=k)):
    kmerdict["".join(kmer)] = float(0)

for line in scores:
    kmer, score = line.rstrip().split("\t")
    kmerdict[kmer] += float(score)

for line in iBED:
    splitline = line.rstrip().split('\t')
    chrom = splitline[0]
    kmer = splitline[3]
    #stt = splitline[2].split('-')[0].split(':')[-1]
    end = splitline[2].split(':')[0]
    stt = splitline[1]
    #end = splitline[2].split('-')[1].split('(')[0]
    #print(splitline)
    #print(stt)
    #print(end)
    #chrom, stt, end, hepta = line.rstrip().split("\t")
    for i in range(0, int(end)-int(stt)):
        if all(c in nct for c in kmer):
            oBED.write("%s\t%s\t%s\t%s\t%.6f\n" % (chrom, str(int(stt)+i), str(int(stt)+i+1), kmer, kmerdict[kmer]))
    #print("")



