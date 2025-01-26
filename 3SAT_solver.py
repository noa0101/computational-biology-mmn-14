'''
This file contains an implementation of the algorithm to solve the SAT problem using molecular computation
'''

from test_tube import Test_Tube
from DNA_strand import DNA_strand
import random
import math

NODE_REP_LENGTH = 20
BASES = ['A', 'T', 'G', 'C']


class SAT_Solver:
    def __init__(self, mnt_of_variables, cnf_representation):
        self.num_of_vars = mnt_of_variables
        self.cnf = cnf_representation
        # literals and nodes get random representations of size NODE_REP_LENGTH
        self.literals_representation = [(''.join(random.choices(BASES, k=NODE_REP_LENGTH)),
                                         ''.join(random.choices(BASES, k=NODE_REP_LENGTH))) for _ in range(mnt_of_variables)]
        self.nodes_representation = [''.join(random.choices(BASES, k=NODE_REP_LENGTH)) for _ in range(mnt_of_variables+1)]

        # edges are represented as described in the algorithm presented in the book
        self.edges_representation = [self.nodes_representation[0]+self.literals_representation[0][0][:NODE_REP_LENGTH//2],
                                     self.nodes_representation[0]+self.literals_representation[0][1][:NODE_REP_LENGTH//2],
                                     self.literals_representation[-1][0][NODE_REP_LENGTH//2:]+self.nodes_representation[-1],
                                     self.literals_representation[-1][1][NODE_REP_LENGTH//2:]+self.nodes_representation[-1]]

        for i in range(1, mnt_of_variables):
            for j in range(2):
                self.edges_representation.append(self.get_edge(self.nodes_representation[i], self.literals_representation[i][j]))

        for i in range(0, mnt_of_variables-1):
            for j in range(2):
                self.edges_representation.append(self.get_edge(self.literals_representation[i][j], self.nodes_representation[i+1]))

        print(
            f'literals representation:\n{self.literals_representation}\n\nnodes_representation:\n{self.nodes_representation}\n')
        print(f'edges representation:\n{self.edges_representation}\n')

        self.test_tube = Test_Tube([])

    # returns the edge between node1 and node2
    @staticmethod
    def get_edge(node1, node2):
        return node1[NODE_REP_LENGTH // 2:] + node2[:NODE_REP_LENGTH // 2]

    # adds to the test tube the representations of the edges
    # and the complementary representations of the literals and nodes
    def prepare_sample(self):
        for edge in self.edges_representation:
            self.test_tube.add_strand(DNA_strand(edge))
        for rep in self.literals_representation:
            for i in range(2):
                self.test_tube.add_strand(DNA_strand(DNA_strand.get_complementary_seq(rep[i])))
        for rep in self.nodes_representation:
            self.test_tube.add_strand(DNA_strand(DNA_strand.get_complementary_seq(rep)))

    # attempts to solve the problem, outputs information from stages in the middle
    def solve(self, pcr_rounds=10):
        print("Attempting to find a solution for the formula: ", self.cnf)
        self.prepare_sample()
        print("size of sample: ", self.test_tube.get_size())
        self.test_tube.pcr(rounds=pcr_rounds)
        print("size first after pcr: ", self.test_tube.get_size())
        self.test_tube.bond()
        print("size after bond: ", self.test_tube.get_size())
        # filter only the bigger strands that have a better chance of becoming entire paths
        self.test_tube.gel_electrophoresis(
            [(self.num_of_vars + 1) * NODE_REP_LENGTH, (2 * self.num_of_vars + 1) * NODE_REP_LENGTH])
        print("size after first gel_electrophoresis: ", self.test_tube.get_size())
        self.prepare_sample()
        self.test_tube.pcr(rounds=pcr_rounds//2)
        print("size after second pcr: ", self.test_tube.get_size())
        self.test_tube.bond()
        print("size after final bonding: ", self.test_tube.get_size())
        # fileter only the strands that represent entire paths in the graph
        self.test_tube.gel_electrophoresis([(2*self.num_of_vars+1)*NODE_REP_LENGTH]*2)
        print("size after final gel_electrophoresis: ", self.test_tube.get_size())
        self.test_tube.pcr(rounds=2)

        # separate only the strands that represent valid solutions
        for C in self.cnf:
            sequences = []
            for x in C:
                sequences.append(self.literals_representation[abs(x)-1][0 if x > 0 else 1])
            self.test_tube.magnetic_separation(sequences)

        print("size after pcr and magnetic separation: ", self.test_tube.get_size())
        solution = self.test_tube.get_one_strand()
        return None if solution is None else self.interpret_solution(solution)

    # takes a DNA strand that represents a valid solution and extracts the corresponding assignment
    def interpret_solution(self, strand):
        assignment = [False]*self.num_of_vars

        # get the correct one of the strands that has the nodes and literals representation
        if strand.sequences[0].startswith(self.nodes_representation[0]):
            seq = strand.sequences[0]
        else:
            seq = DNA_strand.reverse(strand.sequences[1])

        # gets the assignment for each variable
        for i in range(self.num_of_vars):
            assignment[i] = True if seq[(2*i+1)*NODE_REP_LENGTH:].startswith(self.literals_representation[i][0]) else False
        return assignment


if __name__ == "__main__":
    # a few examples
    cnf4 = [[1, 2, 4], [-1, -2, 3], [2, -3, 4], [-1, 2, 4], [2, 3, -4]]
    cnf5 = [[1, 2, 3], [2, -4, 5], [1, -3, 4], [2, 4, -5], [-1, 2, -3]]
    cnf6 = [[1, 3, 5], [2, -4, 6], [1, 3, 6], [-1, 2, -6], [1, 4, -5], [-3, 5, -6]]
    cnf7 = [[1, 4, 7], [5, -6, 7], [2, 3, -5], [-1, -3, -4], [2, -5, 6], [3, 4, -7]]
    cnf8 = [[1, -2, 5], [-3, 4, -6], [2, -4, 7], [1, -5, -8], [-2, 3, -7], [4, -6, 8]]

    solver = SAT_Solver(6, cnf6)
    solution = solver.solve(pcr_rounds=8)
    if solution is None:
        print("No solution found.")
    else:
        print("solution: ",  solution)
