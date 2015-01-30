#!/usr/bin/python
import sys

### This is a throwaway: clean up a reference seed file
### Makes up for the shortcomings in R matrix indexing that Python does not
### have
SEP = '\t'

if __name__ == '__main__':
    with open(sys.argv[1]) as infile:
        line = infile.readline().strip()
        print line
        for line in infile:
            row = line.strip().split(SEP)
            if len(row) > 0:
                row[0] = row[0].replace('"', '').replace('.', '-')
                row[0] = row[0].replace('-txt', '.txt')
                if row[0].startswith("X"):
                    row[0] = row[0].lstrip("X")
                print SEP.join(row)
            
