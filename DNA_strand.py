'''
This file contains the implementation of a single DNA strand, as well as methods that operate locally on given strands
'''

import random

SUCCESS_CHANCE = 0.9  # The success chance of any genetic operation, such as bonding, magnetic separation...


class DNA_strand:
    MATCH = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self, sequence1='', sequence2=None, start=0):
        if DNA_strand.is_legal(sequence1, sequence2, start):
            '''
            The first sequence is always represented from end 5' to end 3', the second in reverse.
            This is taken into account when connecting sequences.
            '''
            self.sequences = [sequence1, sequence2]  # in case of a single strand, sequences[1] is None
            self.start = start  # starting point of the second sequence in comparison to the first

        else:  # DNA strand given is illegal
            raise ValueError('Impossible parameters for DNA strand.')

    # returns whether the strand given describes a legal DNA strand.
    @staticmethod
    def is_legal(sequence1, sequence2, start):
        if sequence2 is None:  # single-stranded
            return sequence1 is not None and len(sequence1) > 0

        # define a minimal number of pairs needed to consider two strands as a double strand
        min_connection = min(len(sequence1), len(sequence2))/2
        if min(len(sequence1), len(sequence2)+start) - max(0, start) < min_connection:
            return False

        for i in range(max(start, 0), min(len(sequence2)+start, len(sequence1))):
            if sequence1[i] != DNA_strand.MATCH[sequence2[i-start]]:  # if pairs don't match
                return False

        return True

    # returns the DNA strand described by parameters of it is legal, otherwise returns None
    @staticmethod
    def get_dna_if_legal(sequence1='', sequence2=None, start=0):
        try:
            return DNA_strand(sequence1, sequence2, start)
        except Exception:
            return None

    # reverses a sequence - useful for converting between 5' to 3' to reverse
    @staticmethod
    def reverse(sequence):
        return ''.join(reversed(sequence))

    # returns the complementary sequence for the given sequence, both represented from 5' to 3' end
    @staticmethod
    def get_complementary_seq(sequence):
        comp = []
        for base in sequence:
            comp.append(DNA_strand.MATCH[base])
        # the complementary strand is from end 3' to end 5'. reverse it to get end 5' to end 3'.
        return DNA_strand.reverse(comp)

    # return the size of the bigger one of the strands composing the DNA piece
    def get_size(self):
        return max(len(self.sequences[0]), (0 if self.sequences[1] is None else len(self.sequences[1])))

    def is_double(self):
        return self.sequences[1] is not None

    # if the DNA strand contains the given sequence, returns the location (strand and index).
    # otherwise returns None, None
    def contains(self, sequence):
        try:
            return 0, self.sequences[0].index(sequence)
        except Exception:
            try:
                return 1, DNA_strand.reverse(self.sequences[1].index(sequence))
            except Exception:
                return None, None

    # performs the action of the described restriction enzyme on the DNA strand
    def restriction_enzyme(self, sequence, dist1, dist2):
        if random.random() > SUCCESS_CHANCE:  # possibly, the action fails
            return
        if self.double:
            strand, index = self.contains(sequence)
            if strand is not None:
                self.sequences[strand] = self.sequences[strand][0:min(index+dist1), len(self.sequences[strand])-1]
                strand ^= 1
                self.sequences[strand] = self.sequences[strand][0:min(index+dist2), len(self.sequences[strand])-1]

    # returns whether the DNA strand contains the given sequence, with a chance of a mistake
    def is_magnetized(self, sequence):
        if random.random() > SUCCESS_CHANCE:  # possibly, the action fails
            return False
        strand, index = self.contains(sequence)
        return strand is not None

    # returns the length difference between the second and first strands at the end
    def get_end(self):
        return (0 if self.sequences[1] is None else len(self.sequences[1])) + self.start - len(self.sequences[0])

    # We will consider two double-stranded DNA's to be compatible for connection iff their sticky ends match exactly
    @staticmethod
    def double_strand_compatability(dna1, dna2):
        if dna2.start == 0 or dna1.get_end() != dna2.start:
            return False
        return DNA_strand.is_legal(dna1.sequences[0]+dna2.sequences[0], dna1.sequences[1]+dna2.sequences[1], dna1.start)

    # returns the result of a connection between the two given strands, or None if it isn't compatible
    @staticmethod
    def connect(dna1, dna2):
        if random.random() > SUCCESS_CHANCE:  # possibly, the action fails
            return None

        possible_dnas = []  # possible connection results (will also contain None for impossible configurations)

        if not (dna1.is_double() or dna2.is_double()):  # two single strands
            seq1 = dna1.sequences[0]
            # reverse one of the strands to ensure connection of end 5' to end 3'
            seq2 = DNA_strand.reverse(dna2.sequences[0])
            for i in range(-1*len(seq1)+1, len(seq2)):
                possible_dnas.append(DNA_strand.get_dna_if_legal(seq1, seq2, i))

        elif dna1.is_double() and dna2.is_double():  # two double strands (possibly can connect at sticky ends)
            # There are four configurations that allow correct ends connection
            if DNA_strand.double_strand_compatability(dna1, dna2):
                possible_dnas.append(DNA_strand(dna1.sequences[0] + dna2.sequences[0], dna1.sequences[1] + dna2.sequences[1],
                                  dna1.start))
            if DNA_strand.double_strand_compatability(dna2, dna1):
                possible_dnas.append(DNA_strand(dna2.sequences[0] + dna1.sequences[0], dna2.sequences[1] + dna1.sequences[1],
                                  dna2.start))

            # flip one dna to explore more connection options
            dna1_flipped = DNA_strand(DNA_strand.reverse(dna1.sequences[1]), DNA_strand.reverse(dna1.sequences[0]),
                                      len(dna1.sequences[1]) + dna1.start - len(dna1.sequences[0]))

            if DNA_strand.double_strand_compatability(dna1_flipped, dna2):
                possible_dnas.append(DNA_strand(dna1_flipped.sequences[0] + dna2.sequences[0],
                                  dna1_flipped.sequences[1] + dna2.sequences[1], dna1_flipped.start))
            if DNA_strand.double_strand_compatability(dna2, dna1_flipped):
                possible_dnas.append(DNA_strand(dna2.sequences[0] + dna1_flipped.sequences[0],
                                  dna2.sequences[1] + dna1_flipped.sequences[1], dna2.start))

        else:  # one double-stranded and the other single-stranded
            if dna2.is_double():
                curr_dna1, curr_dna2 = dna2, dna1
            else:
                curr_dna1, curr_dna2 = dna1, dna2

            # curr_dna1 is double-stranded and curr_dna2 is single-stranded

            if curr_dna1.get_end() > 0:
                possible_dnas.append(DNA_strand.get_dna_if_legal(curr_dna1.sequences[0] + curr_dna2.sequences[0], curr_dna1.sequences[1], curr_dna1.start))

            elif curr_dna1.get_end() < 0:
                possible_dnas.append(DNA_strand.get_dna_if_legal(curr_dna1.sequences[0], curr_dna1.sequences[1] + DNA_strand.reverse(curr_dna2.sequences[0]), curr_dna1.start))

            if curr_dna1.start < 0:
                possible_dnas.append(DNA_strand.get_dna_if_legal(curr_dna2.sequences[0] + curr_dna1.sequences[0], curr_dna1.sequences[1],
                                      len(curr_dna2.sequences[0])+curr_dna1.start))

            elif curr_dna1.start > 0:
                possible_dnas.append(DNA_strand.get_dna_if_legal(curr_dna1.sequences[0], DNA_strand.reverse(curr_dna2.sequences[0]) + curr_dna1.sequences[1],
                                      curr_dna1.start - len(curr_dna2.sequences[0])))

        # if we have obtained a connection possibility (that is not None), return it.
        for dna in possible_dnas:
            if dna is not None:
                return dna

        # otherwise, connection is not possible, return None
        return None

    # these following two methods make sure that a DNA_strand can be nicely printed
    def __str__(self):
        stri = '\n'
        if self.start < 0:
            stri += ' '*abs(self.start)
        stri += self.sequences[0]
        stri += '\n'
        if self.start > 0:
            stri += ' '*self.start
        stri += '' if self.sequences[1] is None else self.sequences[1]
        stri += '\n'
        return stri

    def __repr__(self):
        return self.__str__()
