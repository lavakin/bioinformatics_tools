#!/usr/bin/env python3

from Bio import pairwise2
import sys

def get_editing_distace(seq1:str, seq2:str):
    return pairwise2.align.globalms(seq1, seq2, 0, -1, -1, -1)

"""
for a in pairwise2.align.globalms(seq1, seq2, 0, -1, -1, -1):
    print("editing distance:" + str((-1)*int(a[2])))
    print("final sequences:\n" + a[0] + "\n" + a[1] + "\n")
"""
