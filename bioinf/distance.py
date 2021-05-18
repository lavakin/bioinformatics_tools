#!/usr/bin/env python3

import sys
from Bio import pairwise2


def get_editing_distace(seq1:str, seq2:str):
    """
    :param seq1: sequence one
    :param seq2: sequence two
    :return: editing distance of two sequences along with all alignments with the maximum score
    """
    return pairwise2.align.globalms(seq1, seq2, 0, -1, -1, -1)


class SequencesNotTheSameLength(Exception):
    def __init__(self, message="Sequences does not have the same length"):
        """
        :param message: error message
        """
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message}'


def hamming_distance(seq1, seq2):
    """
    :param seq1: sequence one
    :param seq2: sequence two
    :return: hamming distance of the two sequences, if they are the same length
    """
    if len(seq1) == len(seq2):
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    else:
        raise SequencesNotTheSameLength()


    
print(hamming_distance(sys.argv[1], sys.argv[2]))

"""
for a in pairwise2.align.globalms(seq1, seq2, 0, -1, -1, -1):
    print("editing distance:" + str((-1)*int(a[2])))
    print("final sequences:\n" + a[0] + "\n" + a[1] + "\n")
"""
