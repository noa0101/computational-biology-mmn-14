from test_tube import Test_Tube
from DNA_strand import DNA_strand
import random
import math
NODE_REP_LENGTH = 4
BASES = ['A', 'T', 'G', 'C']


class SAT_Solver:
    def __init__(self, mnt_of_variables, cnf_representation):
        self.num_of_vars = mnt_of_variables
        self.cnf = cnf_representation
        self.literals_representation = [(''.join(random.choices(BASES, k=NODE_REP_LENGTH)),
                                         ''.join(random.choices(BASES, k=NODE_REP_LENGTH))) for _ in range(mnt_of_variables)]
        self.nodes_representation = [''.join(random.choices(BASES, k=NODE_REP_LENGTH)) for _ in range(mnt_of_variables+1)]
        self.edges_representation = [self.nodes_representation[0]+self.literals_representation[0][0][:NODE_REP_LENGTH//2],
                                     self.nodes_representation[0]+self.literals_representation[0][1][:NODE_REP_LENGTH//2],
                                     self.literals_representation[-1][0][NODE_REP_LENGTH//2:]+self.nodes_representation[-1],
                                     self.literals_representation[-1][1][NODE_REP_LENGTH//2:]+self.nodes_representation[-1]]

        if mnt_of_variables > 1:
            for j in range(2):
                self.edges_representation.append(self.get_edge(self.literals_representation[0][j], self.nodes_representation[1]))

        for i in range(1, mnt_of_variables):
            for j in range(2):
                self.edges_representation.append(self.get_edge(self.nodes_representation[i], self.literals_representation[i][j]))

        print(
            f'literals representation:\n{self.literals_representation}\n\nnodes_representation:\n{self.nodes_representation}\n')
        print(f'edges representation:\n{self.edges_representation}\n')

        self.test_tube = Test_Tube([])

    def prepare_sample(self):
        for edge in self.edges_representation:
            self.test_tube.add_strand(DNA_strand(edge))
        for rep in self.literals_representation:
            for i in range(2):
                self.test_tube.add_strand(DNA_strand(DNA_strand.get_complementary_seq(rep[i])))
        for rep in self.nodes_representation:
            self.test_tube.add_strand(DNA_strand(DNA_strand.get_complementary_seq(rep)))

    def solve(self):
        self.prepare_sample()
        print("size before pcr: ", self.test_tube.get_size())
        self.test_tube.pcr(rounds=5*self.num_of_vars)
        print("size after pcr: ", self.test_tube.get_size())
        print(set(self.test_tube.sample))
        self.test_tube.bond()
        print("size after bond: ", self.test_tube.get_size())
        print(self.test_tube.sample)
        self.test_tube.gel_electrophoresis([2*self.num_of_vars+1]*2)
        print("size after gel_electrophoresis: ", self.test_tube.get_size())
        for C in self.cnf:
            sequences = []
            for x in C:
                sequences.append(self.literals_representation[abs(x)-1][0 if x > 0 else 1])
            self.test_tube.magnetic_separation(sequences)

        solution = self.test_tube.get_one_strand()
        return None if solution is None else self.interpret_solution(solution)

    def interpret_solution(self, strand):
        assignment = [False]*self.num_of_vars
        seq = ''
        if strand.sequences[0].startswith(self.nodes_representation[0]):
            seq = strand.sequences[0]
        else:
            seq = DNA_strand.reverse(strand.sequences[1])
        for i in range(self.num_of_vars):
            assignment[i] = True if seq[(2*i+1)*NODE_REP_LENGTH].startswith(self.literals_representation[i][0]) else False
        return assignment

    @staticmethod
    def get_edge(node1, node2):
        return node1[NODE_REP_LENGTH//2:] + node2[:NODE_REP_LENGTH//2]


if __name__ == "__main__":
    cnf = [[1, 2, -3], [-2, 3], [-1, 3]]
    cnf1 = [[1, -1], [1]]
    solver = SAT_Solver(1, cnf1)
    solution = solver.solve()
    print("No solution found" if solution is None else solution)
