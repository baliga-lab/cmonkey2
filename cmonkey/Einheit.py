'''
Created on Mar 26, 2012

@author: frank schmitz, sbri, 2012 -
'''

DNAcD = {"A" : "T", "C" : "G", "G" : "C", "T" : "A", "N" : "N"}


class Einheit():
    """
    the Einheit (unit) represents a gene, gene product, class of genes
    it is the highest unit and corresponds to the matrix rows (eg genes)
    
    
    """
    def __init__(self, Name, Description, GeneID):
        self.__GeneName = Name
        self.__GeneDescription = Description
        self.__GeneID = GeneID
        self.__Subunits = []
        self.__SeqLD = {}
        self.__BSeqLD = {}
     
    def getGeneDescription(self):
        return self.__GeneDescription
    def getGeneID(self):
        return self.__GeneID
    def getGeneName(self):
        return self.__GeneName
    
    def getGeneContig(self):
        return self.__contig
    def setContig(self, contig):
        self.__contig = contig


    def pushSubUnit(self, subUnit):
        self.__Subunits.append(subUnit)
    def getSubUnits(self):
        return self.__Subunits

    def getSeqL(self, keyword):
        try:
            tempL = self.__SeqLD[keyword]
        except KeyError:
            tempL = []
        return tempL
    def getSeqDs(self):
        return (self.__SeqLD, self.__BSeqLD)
    
    def putSeq(self, keyword, Seq):
        try:
            tempL = self.__SeqLD[keyword]
        except KeyError:
            tempL = []
        if Seq not in tempL:
            tempL.append(Seq)
        self.__SeqLD[keyword] = tempL
        return len(tempL)

    def getSeqKeywords(self):
        return self.__SeqLD.keys()



    def updateEinheit_Seq(self, content):
        """
        a huge string is delivered with the entirety of DNA seq
        
        """

        for SubU in self.__Subunits:
            kws = SubU.getSeqKeywords()
            for kw in kws:
                tempSeqL = []
                tempBSeqL = []
                
                SeqL = SubU.getSeqs(kw)
                for Seq in SeqL:
                    Seq.rebirthSeq(content)
                    tempSeqL.append(Seq.getSeq())
                    tempBSeqL.append(Seq.getBSeq())
                
                                                    # push in the sequences per subunit per Unit 
                try:
                    tempL = self.__SeqLD[kw]
                except KeyError:
                    tempL = []
                self.__SeqLD[kw] = tempL + tempSeqL 
                try:
                    tempL = self.__BSeqLD[kw]
                except KeyError:
                    tempL = []
                self.__BSeqLD[kw] = tempL + tempBSeqL
                    





        
class SubEinheit():
    """
    The SubEinheit (SubUnit) represents something like
    as protein isoform or transcripts
    it still can contain multiple sequences of different classes
    
    """
    
    def __init__(self, Name):
        self.__SubID = Name
        self.__species = ""
        self.__SeqD = {}

     
    def getName(self):
        return self.__SubID

    def pushSeq(self, contig, Start, Stop, Strand, SubType, boundaries):
        try:
            tempL = self.__SeqD[SubType]
        except:
            tempL = []
        tempL.append(Sequence(contig, Start, Stop, Strand, SubType, boundaries))
        self.__SeqD[SubType] = tempL
    def replaceSeqL(self, SubType, SeqL):
        self.__SeqD[SubType] = SeqL

    def getSeqs(self, keyword = None):
        if keyword == None:
            return self.__SeqD
        else:
            try:
                return self.__SeqD[keyword]
            except KeyError:
                return None
    def getSeqKeywords(self):
        return self.__SeqD.keys()





class Sequence():
    
    def __init__(self, contig, Start, Stop, Strand, SubType, boundaries):
        self.__contig = contig
        self.__Seq = ""
        self.__BSeq = ""
        self.__Start = Start
        self.__Stop = Stop
        self.__Strand = Strand
        self.__Type = "Generic"
        self.__SubType = SubType
        self.__boundaries = boundaries
        self.__reincarn = "earth"
    def getStart(self):
        return self.__Start
    def getStop(self):
        return self.__Stop
    def getStrand(self):
        return self.__Strand
    def getType(self):
        return self.__Type
    def getSubType(self):
        return self.__SubType
    def getSeq(self):
        return self.__Seq
    def getBSeq(self):
        return self.__BSeq
    def getSeqs(self):
        return (self.__Seq, self.__BSeq)
    def rebirthSeq(self, content):
        if self.__Strand.upper() == "R":
            ref = self.__Stop
            mult = -1
        else:
            ref = self.__Start
            mult = 1
        start = ref + (mult * self.__boundaries[0])
        stop = ref + (mult * self.__boundaries[1])
        bstart = ref + (mult * self.__boundaries[2])
        bstop = ref + (mult * self.__boundaries[3])
        
        bracket = content[start : stop]
        bbracket = content[bstart : bstop]
        if self.__Strand.upper() == "R":
            bracket = "".join([DNAcD[x] for x in reversed(bracket)])
            bbracket = "".join([DNAcD[x] for x in reversed(bbracket)])
        self.__Seq = bracket
        self.__BSeq = bbracket
        self.__reincarn = "Nirvana"
    def getRebirth(self):
        return self.__reincarn
    
    
        
        

