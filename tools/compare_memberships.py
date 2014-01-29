#!/usr/bin/python
import sys

def compare(file1, file2):
    with open(file1) as infile1:
        with open(file2) as infile2:
            header1 = infile1.readline().strip().split('\t')
            header2 = infile2.readline().strip().split('\t')
            if len(header1) != len(header2):
                raise Exception('Num clusters do not match')
            inrows1 = [line.strip().split('\t') for line in infile1]
            inrows2 = [line.strip().split('\t') for line in infile2]
            data1 = {line[0]: line[1:] for line in inrows1 }
            data2 = {line[0]: line[1:] for line in inrows2 }
            if len(data1) != len(data2):
                raise Exception('Numbers of entries does not match')
            if set(data1.keys()) != set(data2.keys()):
                print data1.keys()
                print data2.keys()
                raise Exception('Keys do not match')
            for key in data1:
                if data1[key] != data2[key]:
                    raise Exception("data for key '%s' does not match" % key)

    print "done, everything matches"

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print "usage: ./compare_memberships.py <file1> <file2>"
    else:
        compare(sys.argv[1], sys.argv[2])

