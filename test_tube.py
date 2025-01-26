'''
This file contains an implementation of a test tube containing multiple DNA strands.
'''

import copy
from DNA_strand import DNA_strand
import random
from collections import deque

PCR_CHANCE = 0.95  # chance that pcr will duplicate some strand


class Test_Tube:
    # sample is a list of the DNA strands currently in the test tube.
    def __init__(self, sample):
        self.sample = sample

    def add_strand(self, strand):
        self.sample.append(strand)

    # returns amount of DNA strands in the test tube
    def get_size(self):
        return len(self.sample)

    def pcr(self, rounds=10):
        for _ in range(rounds):
            added_list = []
            for strand in self.sample:
                if random.random() < PCR_CHANCE:
                    added_list.append(copy.deepcopy(strand))
            self.sample += added_list

    # magnetically separates all strands containing al least one of the sequences in the parameter "sequences"
    def magnetic_separation(self, sequences):
        new_sample = []
        for strand in self.sample:
            for sequence in sequences:  # check if the strand contains any of the sequences received
                if strand.is_magnetized(sequence):
                    new_sample.append(strand)
                    break
        self.sample = new_sample

    # separates all the strands whose size is within the given range
    def gel_electrophoresis(self, range):
        new_sample = []
        for strand in self.sample:
            if range[0] <= strand.get_size() <= range[1]:
                new_sample.append(strand)
        self.sample = new_sample

    # performs the action of the restriction enzyme described by the parameters on each strand in the test tube
    def restriction_enzyme(self, sequence, dist1, dist2):
        for strand in self.sample:
            strand.restriction_enzyme(sequence, dist1, dist2)

    # bonds the DNA strands that are compatible, in a random fation.
    # in a "wet" lab this happens naturally
    def bond(self):
        rounds = 10  # how many rounds in which no new bonds are created we should wait before stopping the process
        changes = [len(self.sample)]*rounds
        loop = 0
        while sum(changes) > 0:  # while some bond has been created in the last <rounds> rounds
            change = 0
            new_sample = []
            random.shuffle(self.sample)  # shuffle to simulate molecules moving in the test tube
            i = 0
            while i < len(self.sample):
                if i == len(self.sample)-1:
                    new_sample.append(self.sample[i])
                    break
                connection = DNA_strand.connect(self.sample[i], self.sample[i+1])  # try to connect adjacent strands
                if connection is None:  # strands are not compatible
                    new_sample.append(self.sample[i])
                    i += 1
                else:  # strands are compatible
                    new_sample.append(connection)
                    i += 2
                    change += 1

            changes[loop % rounds] = change  # update changes array
            self.sample = new_sample
            loop += 1

    # returns one strand from the tube if it is not empty, else None
    def get_one_strand(self):
        if len(self.sample) == 0:
            return None
        else:
            return self.sample[0]
