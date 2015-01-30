#!/usr/bin/python
import sys
import math

# This is a quick helper script to preprocess the STRING files that come out of
# cMonkey-R
def identity(param):
        return param

def normalize_edge_list(edges, max_score):
    """normalize scores to 1000, for combined scores"""
    result = []
    for edge in edges:
        score = edge[2] / max_score * 1000.0
        score = 1000 * math.exp(score / 1000.0) / math.exp(1.0)
        result.append((edge[0], edge[1], score))
    return result

def normalize(filename, sub=identity):
        print "processing file '%s'" % filename
        with open(filename) as infile:
                lines = infile.readlines()
        rows = [line.strip().split(' ') for line in lines]
        rows = [(sub(row[0]), sub(row[1]), float(row[2])) for row in rows]
        max_score = 0.0
        for row in rows:
                if row[2] > max_score:
                        max_score = row[2]
        rows = normalize_edge_list(rows, max_score)
        for row in rows:
                print "%s\t%s\t%f" % (row[0], row[1], row[2])

if __name__ == '__main__':
        if len(sys.argv) <= 2:
                print "usage: normalize_string <string-file> <prefix>"
        else:
                normalize(sys.argv[1], lambda x: x.replace(sys.argv[2], ''))
