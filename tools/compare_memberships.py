#!/usr/bin/python
import sys
import argparse
import util


def compare(file1, file2, rnames=False):
    num_errors = 0
    with open(file1) as infile1:
        with open(file2) as infile2:
            header1 = infile1.readline().strip().split('\t')
            header2 = infile2.readline().strip().split('\t')
            if len(header1) != len(header2):
                raise Exception('Num clusters do not match')
            inrows1 = [line.strip().split('\t') for line in infile1]
            inrows2 = [line.strip().split('\t') for line in infile2]
            data1 = {util.make_rname(line[0], rnames): set(line[1:]) for line in inrows1}
            data2 = {util.make_rname(line[0], rnames): set(line[1:]) for line in inrows2}
            if len(data1) != len(data2):
                raise Exception('Numbers of entries does not match')
            if set(data1.keys()) != set(data2.keys()):
                for key1 in data1:
                    if key1 not in data2:
                        print "Key '%s' not found" % key1
                raise Exception('Keys do not match')
            for key in data1:
                if data1[key] != data2[key]:
                    print("data for key '%s' does not match" % key)
                    for d1 in data1[key]:
                        if d1 not in data2[key]:
                            num_errors += 1
    if num_errors == 0:
        print "done, everything matches"
    else:
        print "%d errors found" % num_errors

if __name__ == '__main__':
    description = "compare_memberships.py - compare memberships from 2 tsv files"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--rnames', required=False, action="store_true", default=False)
    parser.add_argument('file1')
    parser.add_argument('file2')
    args = parser.parse_args()
    compare(args.file1, args.file2, args.rnames)

