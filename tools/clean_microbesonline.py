#!/usr/bin/python
import csv
import sys
import os

def read_include_conds(filename):
  with open(filename) as infile:
    conds = [line.rstrip() for line in infile.readlines()]
  return conds

if __name__ == '__main__':
  if len(sys.argv) == 1:
    print 'please provide input file'
  else:
    result = []
    specify_keep_cols = len(sys.argv) > 2
    with open(sys.argv[1]) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
        # cut out the locus id column
        for row in reader:
            new_row = [row[0]]
            new_row.extend(row[2:])
            result.append(new_row)
    
    if specify_keep_cols:
      include_conds = read_include_conds(sys.argv[2])
      keep_indexes = [0]
      keep_indexes.extend([result[0].index(name) for name in include_conds])
      new_result = []
      for row in result:
        new_result.append([row[i] for i in keep_indexes])
      #print "KEEP: ", keep_indexes
      #print new_result
      result = new_result

    # remove duplicate columns
    removecols = set()
    for colnum1 in range(len(result[0])):
        if colnum1 not in removecols:
            values1 = [result[rownum][colnum1] for rownum in range(1, len(result))]
            for colnum2 in range(colnum1 + 1, len(result[0])):
                if colnum2 != colnum1:
                    values2 = [result[rownum][colnum2] for rownum in range(1, len(result))]
                    if values1 == values2:
                        removecols.add(colnum2)

    #print "# columns to remove: %d" % len(removecols)
    #print removecols
    # output the results
    new_result = []
    for row in result:
        newrow = [row[colnum] for colnum in range(len(row))
                  if not colnum in removecols]
        print '\t'.join(newrow)
