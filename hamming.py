import sys

class SeqencesNotTheSameLength(Exception):
    def __init__(self, message="Sequences does not have the same length"):
        self.message = message
        super().__init__(self.message)

    def __str__(self):
        return f'{self.message}'


def hamming_distance(seq1, seq2):
    if len(seq1) == len(seq2):
        return sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    else:
        raise SeqencesNotTheSameLength()

print(hamming_distance(sys.argv[1], sys.argv[2]))
