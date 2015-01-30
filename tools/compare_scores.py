#!/usr/bin/python
import argparse
import math
import toolsutil as util

EPS = 0.0001
def compare(file1, file2, verbose, rnames, mapheaders, eps=EPS):
    def tofloat(x):
        return float('nan') if x == 'NA' else float(x)

    num_errors = 0
    num_correct = 0
    
    with open(file1) as infile1:
        with open(file2) as infile2:
            header1 = infile1.readline().strip().split('\t')
            header2 = infile2.readline().strip().split('\t')
            if len(header1) != len(header2):
                raise Exception('Num clusters do not match')
            inrows1 = [line.strip().split('\t') for line in infile1]
            inrows2 = [line.strip().split('\t') for line in infile2]
            data1 = {util.make_rname(line[0], rnames): map(tofloat, line[1:]) for line in inrows1 }
            data2 = {util.make_rname(line[0], rnames): map(tofloat, line[1:]) for line in inrows2 }
            if len(data1) != len(data2):
                raise Exception('Numbers of entries does not match')
            if set(data1.keys()) != set(data2.keys()):
                print data1.keys()
                print data2.keys()
                raise Exception('Keys do not match')
            for key in data1:
                values1 = data1[key]
                values2 = data2[key]
                #print "VALUES1 = ", values1
                #print "VALUES2 = ", values2
                if len(values1) != len(values2):
                    raise Exception("data for key '%s' does not have the same length" % key)
                for i1 in range(len(header1)):
                    if mapheaders:
                        i2 = header2.index(header1[i1])
                    else:
                        i2 = i1
                    if math.isnan(values1[i1]) and  math.isnan(values2[i2]):
                        continue
                    elif math.isnan(values1[i1]) and not math.isnan(values2[i2]):
                        if verbose:
                            print "[%s, %d]: NaN != %.13f" % (key, i1, values2[i2])
                        num_errors += 1
                    elif not math.isnan(values1[i1]) and math.isnan(values2[i2]):
                        if verbose:
                            print "[%s, %d]: %.13f != NaN" % (key, i1, values1[i1])
                        num_errors += 1
                    elif abs(values1[i1] - values2[i2]) > eps:
                        if verbose:
                            print "[%s, %d/%d]: %.13f != %.13f" % (key, i1, i2,
                                                                   values1[i1], values2[i2])
                        num_errors += 1
                        #raise Exception("key '%s' col %d mismatch (%f != %f)" % (key, i, values1[i], values2[i]))
                    else:
                        num_correct += 1

    return num_correct, num_errors

if __name__ == '__main__':
    description = "compare_scores.py - compare the scores from 2 tsv files"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--verbose', required=False, action="store_true", default=False)
    parser.add_argument('--rnames', required=False, action="store_true", default=False)
    parser.add_argument('--mapheaders', required=False, action="store_true", default=False)
    parser.add_argument('file1')
    parser.add_argument('file2')
    args = parser.parse_args()

    num_correct, num_errors = compare(args.file1, args.file2, args.verbose, args.rnames,
                                      args.mapheaders)
    print "done, %d correct, %d errors" % (num_correct, num_errors)
