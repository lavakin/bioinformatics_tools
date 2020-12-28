from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys

for a in pairwise2.align.globalms(sys.argv[1], sys.argv[2], 0, -1, -1, -1):
    print("editing distance:" + str((-1)*int(a[2])))
    print("final sequences:\n" + a[0] + "\n" + a[1] + "\n")
