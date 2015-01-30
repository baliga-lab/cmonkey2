#!/usr/bin/python
import csv
import sys

if __name__ == '__main__':
  if len(sys.argv) == 1:
      print 'please provide input file'
  else:
      with open(sys.argv[1]) as csvfile:
          line = 0
          reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
          for row in reader:
              if line > 0:
                  try:
                      nums = [float(item) for item in row[1:]]
                  except:
                      print 'Number format error in line %d' % (line + 1)
              line += 1

