#!/usr/bin/env python3

import fastaparser
import sys


class IdNotFoundException(Exception):

    def __init__(self, message="Id not found"):
        """
        :param message: error message
        """
        self.message = message
        super().__init__(self.message)


class Sequence(fastaparser.FastaSequence):
    
    def __len__(self):
        return len(self.sequence_as_string())
    
    
    def __getitem__(self, subscript):
        if isinstance(subscript, slice):
            return self.sequence_as_string()[subscript]
        else:
            return self.sequence_as_string()[subscript]

    
    def print_length(self):
        """
        :return: length of the sequence
        """
        print(self.id, len(self.sequence_as_string()))


    def return_subseq(self, lb: int, ub: int):
        """
        :param lb: start index
        :param ub: stop index + 1
        :return: subsequence of the sequence
        """
        return (self.sequence_as_string()[lb:ub])


class Fasta:
    def __init__(self, fasta_file):
        """
        :param fasta_file: name of the fasta file
        """
        with open(fasta_file) as fasta_file:
            self.sequences = [x for x in fastaparser.Reader(fasta_file)]
            for x in self.sequences:
                x.__class__ = Sequence # if a misterious bug occurs, this is the most probable cause


    def get_ids(self):
        """
        :return: ids of the fasta files
        """
        return [x.id for x in self.sequences]


    def get_sequence_by_index(self, index: int):
        """
        :param index: index of the sequence
        :return: sequence on given index
        """
        return self.get_sequence(self.get_ids()[index])


    def get_sequence(self, name):
        """
        :param name:  name of the sequence
        :return: sequence with the given name
        """
        first = next(filter(lambda x: x.id == name, self.sequences), None)
        if first == None:
            raise IdNotFoundException()
        return first


fasta = Fasta(sys.argv[1])
print((fasta.sequences[0][3]))

