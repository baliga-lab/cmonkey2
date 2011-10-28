"""human.py - cMonkey human specific module

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import util
import numpy
import scipy


def genes_per_group_proportional(num_genes_total, num_per_group):
    """takes the total number of genes and a dictionary containing
    the number of members for each group and returns a dictionary that
    distributes the number of genes proportionally to each group"""
    result = {}
    num_group_elems = sum(num_per_group.values())
    groups = num_per_group.keys()
    for index in range(len(groups)):
        group = groups[index]
        if index == len(groups) - 1:
            result[group] = num_genes_total - sum(result.values())
        else:
            result[group] = int(float(num_genes_total) *
                                float(num_per_group[group]) /
                                float(num_group_elems))
    return result


def genes_per_group_nonproportional(num_genes_total, groups):
    """distributes the number of genes evenly to each group given"""
    result = {}
    partition = int(float(num_genes_total) / float(len(groups)))
    for index in range(len(groups)):
        if index == len(groups) - 1:
            result[groups[index]] = num_genes_total - sum(result.values())
        else:
            result[groups[index]] = partition
    return result


def select_probes(matrix, num_genes_total, column_groups, proportional=True):
    """select probes proportional, column_groups is a map from a group
    label to column indexes in the matrix"""
    def coeff_var(row_values):
        """computes the coefficient of variation"""
        sigma = util.r_stddev(row_values)
        mu = numpy.mean(row_values)
        return sigma / mu

    num_per_group = {group: len(indexes)
                     for group, indexes in column_groups.items()}
    if proportional:
        per_group = genes_per_group_proportional(num_genes_total,
                                                 num_per_group)
    else:
        per_group = genes_per_group_proportional(num_genes_total,
                                                 column_groups.keys())

    cvrows = []
    for group, col_indexes in column_groups.items():
        group_cvs = []
        for row in range(matrix.num_rows()):
            row_values = [matrix[row][col] for col in col_indexes]
            group_cvs.append(coeff_var(row_values))
        cvrows += [group_cvs.index(value)
                   for value in
                   sorted(group_cvs, reverse=True)][:per_group[group]]
    return sorted(list(set(cvrows)))


def intensities_to_ratios(matrix, controls, column_groups):
    """turn intensities into ratios
    Warning: input matrix is modified !!!"""
    control_indexes = [matrix.column_names().index(control)
                       for control in controls
                       if control in matrix.column_names()]
    for group_columns in column_groups.values():
        group_controls = [index for index in control_indexes
                          if index in group_columns]
        means = []
        for row in range(matrix.num_rows()):
            values = [float(matrix[row][col]) for col in group_controls]
            means.append(sum(values) / float(len(values)))

        for col in group_columns:
            for row in range(matrix.num_rows()):
                matrix[row][col] /= means[row]

        center_scale_filter(matrix, group_columns, group_controls)
    return matrix


def center_scale_filter(matrix, group_columns, group_controls):
    """center the values of each row around their median and scale
    by their standard deviation. This is a specialized version"""
    centers = [scipy.median([matrix[row][col]
                             for col in group_controls])
               for row in range(matrix.num_rows())]
    scale_factors = [util.r_stddev([matrix[row][col]
                                    for col in group_columns])
                     for row in range(matrix.num_rows())]
    for row in range(matrix.num_rows()):
        for col in group_columns:
            matrix[row][col] -= centers[row]
            matrix[row][col] /= scale_factors[row]
    return matrix


##################
# Organism interface

class Human:

    def __init__(self, prom_seq_filename, p3utr_seq_filename,
                 nw_factories):
        self.__prom_seq_filename = prom_seq_filename
        self.__p3utr_seq_filename = p3utr_seq_filename
        self.__nw_factories = nw_factories
        self.__networks = None

    def networks(self):
        """returns this organism's networks"""
        if self.__networks == None:
            self.__networks = []
            for make_network in self.__nw_factories:
                self.__networks.append(make_network(self))
        return self.__networks
