#!/usr/bin/python
import sys
import csv

if __name__ == '__main__':
    if len(sys.argv) == 1:
        print 'please provide input file'
    else:
        with open(sys.argv[1]) as csvfile:
            reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
            rows = [row for row in reader]
        print "\t".join(rows[0])
        for i in range(1, len(rows)):
            row = rows[i]
            newrow = [row[0].replace('BT', 'BT_')]
            newrow.extend(row[1:])
            print "\t".join(newrow)
