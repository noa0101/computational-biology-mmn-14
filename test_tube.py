import copy
from DNA_strand import DNA_strand
import random

PCR_CHANCE = 0.5
BONDING_STOP_COND = 0.0001


class Test_Tube:
    # sample is a list of the DNA strands currently in the test tube.
    def __init__(self, sample):
        self.sample = sample

    def add_strand(self, strand):
        self.sample.append(strand)

    def get_size(self):
        return len(self.sample)

    def pcr(self, rounds=10):
        for _ in range(rounds):
            added_list = []
            for strand in self.sample:
                if random.random() < PCR_CHANCE:
                    added_list.append(copy.deepcopy(strand))
            self.sample += added_list

    def magnetic_separation(self, sequences):
        new_sample = []
        for strand in self.sample:
            for sequence in sequences:
                if strand.is_magnetized(sequence):
                    new_sample.append(strand)
                    break
        self.sample = new_sample

    def gel_electrophoresis(self, range):
        new_sample = []
        for strand in self.sample:
            if range[0] <= strand.get_size() <= range[1]:
                new_sample.append(strand)
        self.sample = new_sample

    def restriction_enzyme(self, sequence, dist1, dist2):
        for strand in self.sample:
            strand.restriction_enzyme(sequence, dist1, dist2)

    def bond(self):
        changes = len(self.sample)
        #while changes > BONDING_STOP_COND * len(self.sample):
        for _ in range(100):
            changes = 0
            new_sample = []
            random.shuffle(self.sample)
            i = 0
            while i < len(self.sample):
                if i == len(self.sample)-1:
                    new_sample.append(self.sample[i])
                    break
                connection = DNA_strand.connect(self.sample[i], self.sample[i+1])
                if connection is None:
                    new_sample.append(self.sample[i])
                    i += 1
                else:
                    print(f"connecting: {self.sample[i]}end={self.sample[i].get_end()}\nwith{self.sample[i+1]}end={self.sample[i+1].get_end()}\n")
                    print(f"to create: {connection}")
                    new_sample.append(connection)
                    i += 2
                    changes += 1
            self.sample = new_sample

    def get_one_strand(self):
        if len(self.sample) == 0:
            return None
        else:
            return self.sample[0]
