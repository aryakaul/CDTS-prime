#!/usr/bin/env python

"""

script used to put all allelic frequency on the same position (with order A C G T).
and keep total AN (allele number) and whether the most frequent variant is a singleton or not
and whether the annotated reference allele is the most prevalent or not among the studied population,
and finally which nucleotide is the most prevalent at the position (new reference).

use : input | python reformat_ACAFAN_1linePerPos_addFlags.py > output

input structure : chromosome; start; end; annotated reference; annotated alt; AC; AF; AN;
output structure: chromosome; start; end; annotated reference; AF_A; AF_C; AF_G; AF_T; nct_withMaxAf; sumAltAf; maxIsSgltYesOrNo; AN
where
    AF or Af stand for allelic frequency
    nct stand for nucleotide
    nct_withMaxAf, is the nucleotide with the highest af in the population (which becomes the reference in the studied population)
    sumAltAf, is the combined allelic frequency of the non-reference alleles
    maxIsSgltYesOrNo, indicates whether the alternative allele with the highest allelic frequency is a singleton or not
    AN stands for allele number

"""

import sys

iFreq = sys.stdin
oFreq = sys.stdout
# iFreq is the input file
# oFreq is the output file

Afreq = 0
Cfreq = 0
Gfreq = 0
Tfreq = 0

nctList=['A','C','G','T']

line = iFreq.readline()
pchrom, pstt, pend, pref, palt, pac, paf, pan = line.rstrip().split("\t")

# chrom stands for chromosome
# stt stands for start
# end stands for end
# ref stands for reference
# alt stands for alternative
# ac stands for allele count
# af stands for allele frequency
# an stands for allele number
# p nomenclature in front of any name is used for "previous"

locals()[pref+"freq"] = 1
locals()[palt+"freq"] = float(paf)
locals()[pref+"freq"] -= float(paf)

oFreq.write("%s\t%s\t%s\t%s\t" % (pchrom, pstt, pend, pref))

line = iFreq.readline()
chrom, stt, end, ref, alt, ac, af, an = line.rstrip().split("\t")

while line :
    while stt == pstt :
        locals()[alt+"freq"] = float(af)
        locals()[ref+"freq"] -= float(af)
        pchrom, pstt, pend, pref, palt, pac, paf, pan = chrom, stt, end, ref, alt, ac, af, an
        line = iFreq.readline()
        if line == "":
            break
        chrom, stt, end, ref, alt, ac, af, an = line.rstrip().split("\t")
    A, C, G, T = float(Afreq), float(Cfreq), float(Gfreq), float(Tfreq)
    idxMax = [A,C,G,T].index(max(A,C,G,T))
    maxNct = nctList[idxMax]
    sumAltAf = 1-float(locals()[maxNct+"freq"])
    if round(sumAltAf * float(pan)) == 1:
        sgltFlag="sglt" # sglt stands for singleton
    else:
        sgltFlag="notSglt"
    oFreq.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
            (str(locals()[pref+"freq"]), str(Afreq), str(Cfreq), str(Gfreq), str(Tfreq), maxNct, str(sumAltAf), sgltFlag, pan))
    oFreq.write("%s\t%s\t%s\t%s\t" % (chrom, stt, end, ref))
    pchrom, pstt, pend, pref, palt, pac, paf, pan = chrom, stt, end, ref, alt, ac, af, an
    Afreq = 0
    Cfreq = 0
    Gfreq = 0
    Tfreq = 0
    locals()[pref+"freq"] = 1
    locals()[palt+"freq"] = float(af)
    locals()[pref+"freq"] -= float(af)
    line = iFreq.readline()
    if line == "":
        break
    chrom, stt, end, ref, alt, ac, af, an = line.rstrip().split("\t")

oFreq.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % \
    (str(locals()[pref+"freq"]), str(Afreq), str(Cfreq), str(Gfreq), str(Tfreq), maxNct, str(sumAltAf), sgltFlag, pan))









