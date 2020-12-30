from Bio import AlignIO
import itertools
from Bio.SubsMat import MatrixInfo

class IdNotFoundException(Exception):

    def __init__(self, message="Id not found"):
        self.message = message
        super().__init__(self.message)


def get_sequence(align, index):
    if type(index) is int:
        return align[index]
    elif type(index) is str:
        first =  next(filter(lambda x: x.id == index, align), None)
        if first == None:
            raise IdNotFoundException()
        else:
            return first
    else:
        raise FormatError()


def get_column(align, index):
    return align[:,index]


def get_sum_of_pairs_column(seqences, col, matrix):
    return get_sum_of_pairs(get_column(seqences, col), matrix)


def get_sum_of_pairs(sequences, matrix):
    pairs = list(itertools.combinations(sequences, 2))
    score = 0
    for pair in pairs:
        score += cal_pairwise_score(pair[0], pair[1], matrix)
    return score

    
# function obtained from https://github.com/PrernaDas/Sum-Of-Pairs-Score-Blosum62.git
def cal_pairwise_score(seq1, seq2, matrix):
    """
    This function constructs pairs taking one element from seq1, and another element from the corresponding position in the second sequence.
    It then searches for the pair in the blosum62 matrix
    and assigns the value/score associted with that match pair to a score variable
    """
    score = 0
    for i in range(len(seq1)):
        if seq1[i] != '-' and seq2[i] != '-': # disregarding pairs in which either or one of the element has '-'
            pair = (seq1[i], seq2[i])
            if pair not in blosum: # the blosum pairs are arranged as ('B', 'A'), ('B', 'C'), so in order to look up for the value of ('A','B'), we need to reverse it ('B', 'A) 
                reverse_pair = tuple(reversed(pair))
                score += blosum[reverse_pair] # score = score + blosum[reverse_pair]
            else:
                score += blosum[pair] # score = score + blosum[pair]
    return score # gives the cumulative score for all the pairs evaluated in range of len(seq1), i.e i takes values from 0 uptil len(seq1)
    
align = AlignIO.read("clust", "clustal")
blosum = MatrixInfo.blosum62
print(get_sum_of_pairs_column(align, 100, blosum))
