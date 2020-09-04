#!/usr/bin/env python

import sys
from statistics import mean

iTXT = sys.stdin
oTXT = sys.stdout
stored_entry = ("NULL", "NULL", "NULL")
cdts_gender_scores = [0,0]
for line in iTXT:
    entry = line.rstrip().split('\t')
    chrom = entry[0]
    st = entry[1]
    en = entry[2]
    cdts = float(entry[3])
    gender = entry[4]
    curr_entry = (chrom, st, en)
    if gender == "MALE": cdts_gender_scores[1] = cdts
    elif gender == "FEMALE": cdts_gender_scores[0] = cdts
    if stored_entry[0] == "NULL":
        stored_entry = curr_entry
        first = True
    if curr_entry == stored_entry and not first:
        if gender == "MALE": cdts_gender_scores[1] = cdts
        elif gender == "FEMALE": cdts_gender_scores[0] = cdts
        #print(cdts_gender_scores)
        oTXT.write("%s\t%s\t%s\t%s\n" % (stored_entry[0], stored_entry[1], stored_entry[2], cdts_gender_scores[1]-cdts_gender_scores[0]))
        cdts_gender_scores = [0,0]
    stored_entry = curr_entry
    first = False
