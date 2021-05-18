#!/usr/bin/env python3

from Bio import AlignIO, Align
import itertools
import math
import numpy as np


class IdNotFoundException(Exception):

    def __init__(self, message="Id not found"):
        """
        :param message: error message
        """
        self.message = message
        super().__init__(self.message)


class MSA:
    def __init__(self, file_name):
        """
        :param file_name: name of the clustal file
        """
        self.alignments = AlignIO.read(file_name, "clustal")


    def get_sequence(self, index):
        """
        :param index: index of a sequence, string or int
        :return: sequence
        """
        if type(index) is int:
            return self.alignments[index]
        elif type(index) is str:
            first = next(filter(lambda x: x.id == index, self.alignments), None)
            return first    
        else:
            raise FormatError()


    def get_column(self, index):
        return self.alignments[:,index]


    def get_conserv_for_seq(self, seq_index, matrix):
        """
        :param seq_index: index of a particular sequence
        :param matrix: scoring matrix
        :return: conservation score for a sequence based on SoP score
        """
        score = 0
        for i in range(len(self.alignments)):
            if i != seq_index:
                score += MSA.cal_pairwise_score(self.alignments[seq_index], self.alignments[i], matrix)
        return score


    def get_N_best(self, N, matrix):
        """
        :param N: number of best columns
        :param matrix: scoring matrix
        :return: N column indexes with highest SoP score together with the respective score
        """
        N = [(0, -math.inf) for i in range(N)]
        for i in range(1,len(self.alignments[0])):
            if i % 100 == 0:
                print(i)
            score = MSA.get_sum_of_pairs_c(self.get_column(i), matrix)
            if score > N[0][1]:
                N[0] = (i, score)
                N.sort(key=lambda x: x[1])
        return N

    
    def get_sum_of_pairs(self, matrix):
        score = 0
        for i in range(len(self.alignments[0])):
            score += MSA.get_sum_of_pairs_c(self.get_column(i), matrix)
        return score
    

    def get_sum_of_pairs_column(self, col, matrix):
        """
        :param col: index of a column
        :param matrix: scoring matrix
        :return: sum of pairs score for a selected column
        """
        return MSA.get_sum_of_pairs_c(self.get_column(col), matrix)


    def get_sum_of_pairs_c(column, matrix):
        """
        :param matrix: scoring matrix
        oaram column: an array of values from a specific column
        :return: sum of pairs for a column
        """
        score = 0
        for pair in itertools.combinations(column, 2):
            score += MSA.cal_pairwise_score(pair[0], pair[1], matrix)
        return score


    def get_N_best_for_sequence(self, N, seq_index, matrix):
        """
        :param N: number of best columns
        :param seq_index: index of the sequence
        :param matrix: scoring matrix
        :return: N columns with highest SoP score
        """
        N = [(0, -math.inf) for i in range(N)]

        for i in range(len(self.alignments[0])):
            col = self.get_column(i)
            score = 0
            for j in range(len(col)):
                if j != seq_index:
                    score += MSA.cal_pairwise_score(col[j], col[seq_index], matrix)
            if score > N[0][1]:
                N[0] = (i, score)
                N.sort(key=lambda x: x[1])
        return N


    # function obtained from https://github.com/PrernaDas/Sum-Of-Pairs-Score-Blosum62.git
    def cal_pairwise_score(seq1, seq2, matrix):
        """
        This function constructs pairs taking one element from seq1, and another element from the corresponding position in the second sequence.
        It then searches for the pair in the blosum62 matrix
        and assigns the value/score associted with that match pair to a score variable
        """
        score = 0
        for i in range(len(seq1)):
            if seq1[i] != '-' and seq2[i] != '-':  # disregarding pairs in which either or one of the element has '-'
                pair = (seq1[i], seq2[i])
                if pair not in matrix:  # the blosum pairs are arranged as ('B', 'A'), ('B', 'C'), so in order to look up for the value of ('A','B'), we need to reverse it ('B', 'A)
                    reverse_pair = tuple(reversed(pair))
                    score += matrix[reverse_pair]  # score = score + blosum[reverse_pair]
                else:
                    score += matrix[pair]  # score = score + blosum[pair]
        return score  # gives the cumulative score for all the pairs evaluated in range of len(seq1), i.e i takes values from 0 uptil len(seq1)


