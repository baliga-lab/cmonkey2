"""extract_ratios.py - extract ratios from a tab separated file,
and converting gene names to VNG names to generate a test data set.
Currently Halobacterium only.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys

def read_oligo_map(filename):
  """Reads an SBEAMS oligo map file"""
  with open(filename) as infile:
    lines = infile.readlines()
    result = {}
    for i in range(1, len(lines)):
      row = lines[i].strip().split('\t')
      result[row[4]] = row[5]
    return result


if __name__ == '__main__':
  if len(sys.argv) <= 3:
    print "usage: python %s <matrix_output> <oligomap> <num columns>"
  else:
    num_cols = int(sys.argv[3])
    oligo_map = read_oligo_map(sys.argv[2])
    with open(sys.argv[1]) as infile:
      lines = infile.readlines()
      header = lines[1].strip().split('\t')
      first_row = header[0]
      for col_index in range(num_cols):
        first_row += ('\tCond%d' % col_index)
      print first_row
      for index in range(2, len(lines) - 1):
        row = lines[index].strip().split('\t')
        out_row = oligo_map[row[0]]
        for col_index in range(num_cols):
          out_row += ('\t%s' % row[2 + col_index])
        print out_row
