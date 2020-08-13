#!/usr/bin/env python

"""

script used to extract heptamers counts from a nucleotide sequence

use : input | python countkmer.py [columnNb] [k] > output

where columnNb indicates in which column (1 based) the sequence is

"""

import sys, itertools

iTXT = sys.stdin
oTXT = sys.stdout
colNb = int(sys.argv[1]) - 1
k = int(sys.argv[2])

kmer_dict={}  # creates a dictionnary that will contain all possible heptamers
nct='ACTG'   # nct stands for nucleotide

for barcode in list(itertools.product(nct, repeat=k)):
    kmer_dict["".join(barcode)] = 0

for line in iTXT :
    seq = line.rstrip().split("\t")[colNb]
    for i in range(0, len(seq)-(k-1)):
        if all(c in nct for c in seq[i:(i+k)]):
            kmer_dict[seq[i:(i+k)]] += 1

for key in kmer_dict:
    oTXT.write("%s\t%s\n" % (key, str(kmer_dict[key])))


