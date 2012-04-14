"""human_test.py - test classes for human module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import unittest
import human
import datamatrix as dm


class HumanTest(unittest.TestCase):  # pylint: disable-msg=R0904
    """Test class for human module"""

    def test_genes_per_group_proportional_one(self):
        """tests the genes_per_group_proportional function"""
        genes_per_group = human.genes_per_group_proportional(2000, {1:1})
        self.assertEquals(2000, genes_per_group[1])

    def test_genes_per_group_proportional_fifty_fifty(self):
        """tests the genes_per_group_proportional function"""
        genes_per_group = human.genes_per_group_proportional(2000, {1:1, 2:1})
        self.assertEquals(1000, genes_per_group[1])
        self.assertEquals(1000, genes_per_group[2])

    def test_genes_per_group_proportional(self):
        """tests the genes_per_group_proportional function"""
        genes_per_group = human.genes_per_group_proportional(2000, {1:2, 2:1})
        self.assertEquals(1333, genes_per_group[1])
        self.assertEquals(667, genes_per_group[2])

    def test_genes_per_group_nonproportional(self):
        """tests the genes_per_group_proportional function"""
        genes_per_group = human.genes_per_group_nonproportional(2000, [1, 2, 3])
        self.assertEquals(666, genes_per_group[1])
        self.assertEquals(666, genes_per_group[2])
        self.assertEquals(668, genes_per_group[3])

    def test_genes_per_group_nonproportional_one(self):
        """tests the genes_per_group_proportional function"""
        genes_per_group = human.genes_per_group_nonproportional(2000, [1])
        self.assertEquals(2000, genes_per_group[1])

    def test_genes_per_group_nonproportional_6(self):
        """tests the genes_per_group_proportional function"""
        genes_per_group = human.genes_per_group_nonproportional(2000, [1, 2, 3, 4, 5, 6])
        self.assertEquals(333, genes_per_group[1])
        self.assertEquals(333, genes_per_group[2])
        self.assertEquals(333, genes_per_group[3])
        self.assertEquals(333, genes_per_group[4])
        self.assertEquals(333, genes_per_group[5])
        self.assertEquals(335, genes_per_group[6])

    def test_select_probes_proportional(self):
        """tests the select_probes_proportional() function"""
        matrix = dm.DataMatrix(2, 2, ['R1', 'R2'], ['C1', 'C2'], [[1, 3], [2, 4]])
        selgenes = human.select_probes(matrix, 10, {1:[0, 1]})
        self.assertEquals([0, 1], selgenes)

    def test_select_probes_proportional(self):
        """tests the select_probes_proportional() function"""
        matrix = dm.DataMatrix(10, 4,
                               [('R%d' % row) for row in range(1, 11)],
                               [('C%d' % col) for col in range(1, 5)],
                               [[1, 11, 21, 31],
                                [2, 12, 22, 32],
                                [3, 13, 23, 33],
                                [4, 14, 24, 34],
                                [5, 15, 25, 35],
                                [6, 16, 26, 36],
                                [7, 17, 27, 37],
                                [8, 18, 28, 38],
                                [9, 19, 29, 39],
                                [10, 20, 30, 40]])
        selgenes = human.select_probes(matrix, 10, {1:[0, 1], 2: [2, 3]})
        self.assertEquals([0, 1, 2, 3, 4], selgenes)

    def test_intensities_to_ratios(self):
        """tests the intensities_to_ratios() function"""
        matrix = dm.DataMatrix(10, 4,
                               [('R%d' % row) for row in range(1, 11)],
                               [('C%d' % col) for col in range(1, 5)],
                               [[1, 11, 21, 31],
                                [2, 12, 22, 32],
                                [3, 13, 23, 33],
                                [4, 14, 24, 34],
                                [5, 15, 25, 35],
                                [6, 16, 26, 36],
                                [7, 17, 27, 37],
                                [8, 18, 28, 38],
                                [9, 19, 29, 39],
                                [10, 20, 30, 40]])
        ratios = human.intensities_to_ratios(matrix, ['C2', 'C4'], {1:range(4)})
        for row in range(ratios.num_rows()):
            self.assertAlmostEquals(-1.549193, ratios[row][0], places=5)
            self.assertAlmostEquals(-0.774597, ratios[row][1], places=5)
            self.assertAlmostEquals(0.0, ratios[row][2], places=5)
            self.assertAlmostEquals(0.774597, ratios[row][3], places=5)
