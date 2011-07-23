import sys

def read_oligo_map(filename):
  with open(filename) as infile:
    lines = infile.readlines()
    result = {}
    for i in range(1, len(lines)):
      row = lines[i].strip().split('\t')
      result[row[4]] = row[5]
    return result


if __name__ == '__main__':
  if len(sys.argv) <= 1:
    print "usage: python %s <matrix_output>"
  else:
    oligo_map = read_oligo_map('halo_oligo.map')
    with open(sys.argv[1]) as infile:
      lines = infile.readlines()
      header = lines[1].strip().split('\t')
      print "%s\tCond1\tCond2" % (header[0])
      for index in range(2, len(lines) - 1):
        row = lines[index].strip().split('\t')
        print "%s\t%s\t%s" % (oligo_map[row[0]], row[3], row[4])
