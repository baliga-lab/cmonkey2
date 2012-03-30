'''
Helper/Prep scripts for FS's cMonkey Python

@author: frank
'''

# FS specific imports

from tools.FileTools import OpenReadFile, OpenWriteFile                 #provide a file(read) handler
from Environment.Wrappers.FileSystem import AutoInitFileStructure       # init my FileStructure


def ReduceSTRING_v1(file, outfile, oldspeciesL):
    """ parse a STRING file and reduce it to the relevant species links
    for now, replace scores with 1000 """
    speciesL = []
    headerscoreline = "score" 
    for sp in oldspeciesL:
        speciesL.append(str(sp))
        # must because the search within string is comparing strings
        # just in case the speclist has numbers
    
    handle = OpenReadFile(file)
    outhandle = OpenWriteFile(outfile)
    header = handle.readline()
    header = header.rstrip()
    newheader = "protein1 protein2 %s\n" % (headerscoreline)
    outhandle.write(newheader)
    for line in handle:
        line  = line.rstrip()
        lines = line.split(" ")
        inN = lines[0].split(".")
        outN = lines[1].split(".")
        scores = line[2:]
        if inN[0] in speciesL or outN[0] in speciesL:
            # found the species at least ONCE
            scoreline = CalcScore(scores)
            newline = " ".join([inN[1], outN[1], scoreline]) + "\n"
            outhandle.write(newline)
    handle.close()
    outhandle.close()
    return
        
def CalcScore(scoreL):
    """ takes a string List of scores as input
        converts is to another string (sep'ed by spaces, its string, eh)
    """
    scoreS = "1000"
    return scoreS


if __name__ == '__main__':
    Dirs = AutoInitFileStructure()
    outfilename = "HumanProteinLinks.txt"
    SpeciesL = [9606]
    ReduceSTRING_v1(Dirs.StringProtLinks, "".join([Dirs.StringPath, outfilename]), SpeciesL)

