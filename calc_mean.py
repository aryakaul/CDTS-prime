#!/usr/bin/env python

import sys
from statistics import mean

iTXT = sys.stdin
oTXT = sys.stdout
prev_entry = ("NULL", "NULL", "NULL")
running_list = []
for line in iTXT:
    entry = line.rstrip().split('\t')
    chrom = entry[0]
    st = entry[1]
    en = entry[2]
    curr_entry = (chrom, st, en)
    if prev_entry[0] == "NULL":
        prev_entry = curr_entry
    if curr_entry == prev_entry:
        running_list.append(float(entry[5]))
    else:
        oTXT.write("%s\t%s\t%s\t%s\n" % (prev_entry[0], prev_entry[1], prev_entry[2], mean(running_list)))
        running_list = [float(entry[5])]
        prev_entry = curr_entry
oTXT.write("%s\t%s\t%s\t%s\n" % (prev_entry[0], prev_entry[1], prev_entry[2], mean(running_list)))
