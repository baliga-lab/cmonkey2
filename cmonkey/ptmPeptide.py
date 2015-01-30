# vi: sw=4 ts=4 et:
"""ptmPeptide.py - processing module for feeding phosphopetides into cMonkey
This should load up peptide expressions and construct the network motifs.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""

import datamatrix
import util
import os
import re
import numpy as np
import network
import collections

class PTMDataFactory:
    """
    An object for reading in the PTM files and converting them into the correct 
    forms for cMonkey.
    """    
    def __init__(self, ptmExprFile='testdata/peptideExpr.augWithExp23.032312.csv'):
        """  Load the data files
        """
        #Read in the PTM file
        self.peptides = self.readExprFile(ptmExprFile)
        
        #Create a standard expression matrix file for 'datamatrix'
        matrixInfo = self.peptides.getExprMatrix()
        nrows = matrixInfo['matrix'].shape[0]
        ncols = matrixInfo['matrix'].shape[1]
        row_names = matrixInfo['rowNames']
        col_names = matrixInfo['colNames']
        
        self.curDM = datamatrix.DataMatrix(nrows=nrows, ncols=ncols, row_names=row_names, col_names=col_names, values=matrixInfo['matrix'])
    
        #Build a separate network for each kinase
        print('Building the network')

        #Load PTM factors
        self.proteinPeptideLUT = self.proteinToPeptidesLUT()  #Requires peptides to be loaded
        ptmFactors = PTMFactors()
        self.networks = ptmFactors.returnNetworks(score = 1, protLUT=self.proteinPeptideLUT)
        #This isn't quite right, next convert from proteins to peptides.

        #Load the motif regulations
        motifFactors = PTMMotifs()
        motifFactors.returnNetworks(score = 0.5, peptides=row_names)
        #NOTE MUST CONVERT PTMMotif protein names to Systematic.  Trouble, as this is yeast-specific
        #  So can I avoid it and just use peptide match here?  Yes I can!  Although network names will be standard
        print('Stopped working here.  Next thing to do, build a returnNetworks for motifFactors.  Then join them')
        
        #Build the motif network -- 
    def getPeptides(self):
        """
        Return the list of peptides
        """
        return (self.peptides)
    
    def proteinToPeptidesLUT(self):
        """
        Build a LUT with proteins as key and peptide list as values
        """
        protPeptideLUT = {}
        pepList = self.getPeptides()
        for peptide in pepList.itervalues():
            for pepprotein in peptide.getProteins():
                if (protPeptideLUT.has_key(pepprotein)):
                    protPeptideLUT[pepprotein].append(peptide.getSequence())
                else:
                    protPeptideLUT[pepprotein] = [peptide.getSequence()]
            
        return(protPeptideLUT)
        
    def returnDM(self):
        """
        Return the data matrix
        """
        return(self.curDM)
        
    def returnNetworks(self):
        """
        Return the networks
        """
        return(self.networks)
        
    def readExprFile(self, fileName):
        """
        Read in an expression list file with the following columns:
            protein - the protein name
            peptide - the phospho peptide (PTMs in lowercase, Amino Acids in upper
            *.exp.* - the expression ratios at various time points
        All other time points will be ignored.
        """
        infile = util.DelimitedFile.read(fileName, has_header=True, quote='\"')
        headerList = infile.header()[0].split(',')
        proteinIdx = headerList.index('protein')
        peptideIdx = headerList.index('peptide')
        
        #Get the first set of expressions that contain with '.exp.'
        exprIdxs = []
        for i in range(len(headerList)):
            if re.search('\.exp\.', headerList[i]):
                    exprIdxs = exprIdxs + [i]
        
        #Load the data into Peptide objects
        peptideList = Peptides()
        for curLine in infile.lines():
            curLine = curLine[0].split(',')
            peptideStr = curLine[peptideIdx] 
            if (not peptideStr in peptideList.keys()):
                curPeptide = Peptide(peptideStr)
                curPeptide.addProtein(curLine[proteinIdx])
                for idx in exprIdxs:
                    curPeptide.addExpr(headerList[idx], curLine[idx])
                peptideList[peptideStr] = curPeptide
            else:
                peptideList[peptideStr].addProtein(curLine[proteinIdx])
        
        #Add in sequence motif information
        
        return(peptideList)


class Peptide:
    """
    A peptide object that can return the sequence of a peptide, the 
    post-translational modification, and the proteins associated with the peptide.
    """
    def __init__(self, modSequence):
        """create a Peptide instance"""
        self._modSequence_ = modSequence #Contains the modification ex. SpNKITITNDKGR
        self._sequence_ = ''.join(e for e in modSequence if e.isupper()) #Sequence without modifications ex. SNKITITNDKGR
        self._proteins_ = [] #The proteins that contain this sequence
        self._expression_ = {} #Bad for memory footprint, good for OO design
        
    def addProtein(self, protein):
        """Add a protein to the list associated with this peptide"""
        self._proteins_ = self._proteins_ + [protein]
        
    def addExpr(self, exprName, expr):
        """Add an expression to the expression hash table"""
        self._expression_[exprName] = expr 
        
    def getModSequence(self):
        return(self._modSequence_)
    
    def getSequence(self):
        return(self._sequence_ )
    
    def getProteins(self):
        return(self._proteins_ )
    
    def getExpression(self):
        return(self._expression_ )
    
    def getExprVect(self):
        """Return the expression values in a list with a numpy.array (vector) at 
        1 for concatenation"""
        exprNames = []
        exprVals = []
        for key in sorted(self._expression_.keys()):
            exprNames = exprNames + [key]
            try:
                newExpr = float(self._expression_[key])
            except Exception:
                newExpr = float('NaN')
            exprVals = exprVals + [newExpr]
        exprVals = np.array(exprVals, dtype=np.float64) #numpy.NaN?        
        return({'exprNames': exprNames, 'exprVals': exprVals})
    
    def isInProtein(self, protein):
        """Test to see if this peptide occurs in a protein"""     
        for curProt in self._proteins_:
            if (curProt == protein):
                return(True)
        return(False)        

class Peptides(dict):
    """
    A dictionary of peptide objects.  Should act like a hash table containing peptide objects.
    Use this and not just a hash table to return a matrix of values and a formatted network?
    But how to check to make sure it only gets 'Peptide' objects?  Overload the add function?
    """
    
    def getExprMatrix(self):
        
        colNames = []
        rowNames = []
        for key in self.keys():
            rowNames = rowNames + [key]
            exprVect  = self[key].getExprVect()
            
            if (len(colNames) == 0):  #First time through
                colNames = exprVect['exprNames']
                matrix = exprVect['exprVals']
            else:
                #Error checking
                if (len(colNames) != len(exprVect['exprNames'])):
                    #pass_ec=False
                    raise Exception('Peptide expression names do not match length')
                else:
                    for i in range(len(colNames)):
                        if (not colNames[i] == exprVect['exprNames'][i]):
                            raise Exception('Peptide expression names do not match names')
                        
                matrix = np.vstack((matrix, exprVect['exprVals']))
        
        return({'matrix': matrix, 'colNames': colNames, 'rowNames': rowNames})

class PTMMotif:
    """
    A PTM factor and its motif.  For example will associate the kinase PKB 
    with the peptide DCEKLRRRFSSLHFMVEVKGD modifying the 11th residue.
    By default (i.e. < 0) the position will pick the middle residue of peptide
    """
    def __init__(self, modifier, peptide, weight=1, position=-1):
        """
        create a Peptide instance
        modifier:  The factor
        peptide:  The target peptide
        weight:  Should be 1 for experimental or a lower number (ex. 0.5) for predicted
        """
        self._modifier_ = modifier #The phosphotase / kinase / etc.
        self._peptide_ = peptide #Sequence without modifications ex. SNKITITNDKGR
        if position > 0:
            self._position_ = position #The position in the peptide of the modified residue
        else:
            self._position_ = int((len(peptide)) / 2)
        self._weight_ = weight
            
    def matchPeptide(self, peptide):
        """
        Return True / False indicating if peptide matches to this motif
        """
        return(self.matchPeptides(peptide=peptide, motif=self._peptide_, motifModLoc=self.position))

    @staticmethod   
    def matchPeptides(peptide="YYFSVFYYEILNSpPEKACSp", motif='FSVFYYEILNSPEKACSLAKT', peptideModLoc=-1, motifModLoc=-1, minOverlap=7):
        """
        Peptide should be a sequence of residues with the modified one at 
        peptideModLoc, positive base 0 indexing OR preferably indicated by lowercase letters
        Motif should be an odd numbered string of residues with the modified one in the middle
        Return True / False indicating if the two peptides match
        Sparsely tested as of 4/3/12
        """
        rv = False
        if not isinstance(peptideModLoc, list):
            peptideModLoc = [peptideModLoc]
        
        if peptideModLoc[0] < 0:
            peptideModLoc = []
            for i in range(len(peptide)):
                if peptide[i].islower():
                    peptideModLoc.append(i - len(peptideModLoc) -1)  #- len(peptideModLoc) removes space for additional modification annotation
        if (motifModLoc < 0):
            motifModLoc=int((len(motif)) / 2)
        
        for pepModLoc in peptideModLoc:
            curpeptide = ''.join(e for e in peptide if e.isupper())
            curmotif = motif
            #Chop to align the start position of the strings
            if pepModLoc > motifModLoc:
                curpeptide = curpeptide[pepModLoc-motifModLoc:]
            else:
                curmotif = curmotif[motifModLoc-pepModLoc:]
                
            #chop to align the end positions
            if len(curpeptide) > len(curmotif):
                curpeptide = curpeptide[:len(curmotif)-len(curpeptide)]
            else:
                curmotif = curmotif[:len(curpeptide)-len(curmotif)]
            
            if (len(curmotif) > minOverlap) and (curpeptide == curmotif):
                rv = True    
            
        return(rv)

class PTMMotifs:
    """
    Load up the PTMMotif files 
    """
    def __init__(self, experimentalFiles=['testdata/PhosphoELM.tab','testdata/Swiss_phos_exp.tab'], predictedFiles=['testdata/Swiss_phos_putative.tab']):
        """
        Constructor
        """
        self.motifDictionary = {}  #Store a list of factor / motif (PTMMotif) objects
        
        for expFile in (experimentalFiles+predictedFiles):
            infile = util.DelimitedFile.read(expFile, has_header=True, quote='\"')
            weight = 0.25 if predictedFiles.count(expFile) > 0 else 0.5 #predicted weight = .25, real = .5
            for curLine in infile.lines():
                if (len(curLine) > 4):  #Sometimes the middle line is split
                    modLine = ''.join([curstr for curstr in curLine[2:-1]])
                    curLine = curLine[0:2] + [modLine] + [curLine[-1]]  #The 0:2 irks me. starting with 0 means second num is number of elements.  BLEH
                modifier=curLine[2]
                
                if modifier.find('(by') == -1: 
                    modifiers = [modifier]
                else:  #There will be a list of modifiers instead
                    modifiers = modifier.partition('(by')[2]
                    modifiers = modifiers[0:modifiers.find(')')].replace('and','').strip().split(' ')
                
                for modifier in modifiers:
                    if (self.motifDictionary.has_key(modifier)):
                        self.motifDictionary[modifier] = self.motifDictionary[modifier] + [PTMMotif(modifier=modifier, peptide=curLine[-1], weight=weight)]
                    else:
                        self.motifDictionary[modifier] = [PTMMotif(modifier=modifier, peptide=curLine[-1], weight=weight)]
            
        #for expFile in predictedFiles:
        #    infile = util.DelimitedFile.read(expFile, has_header=True, quote='\"')
        #    for curLine in infile.lines():
        #        key = curLine[-1]  #This will simply record the peptide
        #        self.motifDictionary[key] = PTMMotif(modifier=curLine[2], peptide=curLine[-1], weight=0.25)
    
    def returnDict(self):
        return(self.motifDictionary)
    
    def returnNetworks(self, peptides, score = 0.5):
        """
        Return the networks
        TOO SLOW!!!  See notebook entry on page 156
        """
        networks = []
        for motifKey in self.motifDictionary.keys():
            factor=motifKey
            peptideDict = {}
            for peptideObj in self.motifDictionary[motifKey]:
                for peptide in peptides:
                    rv = PTMMotif.matchPeptides(peptide=peptide, motif=peptideObj._peptide_, peptideModLoc=-1, motifModLoc=-1, minOverlap=7)
                    if (rv == True):
                        peptideDict[peptide] = peptideObj._weight_
            #Store the edges
            edgeList = []
            if (len(peptideDict) > 1):
                keys = peptideDict.keys()
                for i in range(len(keys)):
                    weight = peptideDict[keys[i]]
                    for j in range(i+1,len(keys)):
                        newEdge = network.NetworkEdge(keys[i], keys[j], score=weight)
                        edgeList.append(newEdge)
        
            if len(edgeList) > 0:
                networks.append(network.Network(name=factor, edges=edgeList, weight=score))
        return(networks)


class PTMFactors:
    """
    Load up the PTMFactor files 
    """
    def __init__(self, ptmregFiles = ['testdata/PhosSubstrate.csv']):
        """
        Constructor
        """
        self.ptmDictionary = {}
            
        for expFile in ptmregFiles:
            infile = util.DelimitedFile.read(filepath=expFile, has_header=True, sep=',')
            for curLine in infile.lines():
                key = curLine[0] + "." + curLine[1] + "." + curLine[2] #This will simply record the peptide
                self.ptmDictionary[key] = PTMFactor(factor=curLine[0], substrate=curLine[1], coregulator=curLine[2])
    
    def returnDict(self):
        return(self.ptmDictionary)
    
    def getFactorsAndSubstrates(self, incCofactor=True):
        """
        Get the list of unique factors stored in the PTMFactor list with their substrates
        NOTE: As of 4/4/12 ignores scores
        """
        factors = {}
        for curPTM in self.ptmDictionary.values():
            key = curPTM._factor_
            if (incCofactor and not(curPTM._coregulator_ == 'NA')):
                key = key + '.' + curPTM._coregulator_
                
            if ( factors.has_key(key) ):
                factors[key] = factors[key] + [curPTM._substrate_]
            else:
                factors[key] = [curPTM._substrate_]
        
        return(factors)
    
    def returnNetworks(self, protLUT, score = 1):
        """
        Convert the list of PTMFactors into a network list with one network per factor
        protLUT is used to convert the substrate into the peptides
        NOTE: As of 4/4/12 ignores scores
        """
        factorDict = self.getFactorsAndSubstrates()
        networks = []
        for factor in factorDict.keys():
            substrates = factorDict[factor]
            edgeList = {}#NetworkEdges, Dictionary to remove repeats
            for i in range(len(substrates)):
                substrate1 = substrates[i]
                if (protLUT.has_key(substrate1)):
                    for peptide1 in protLUT[substrate1]:
                        for j in range(i+1, len(substrates)):
                            substrate2 = substrates[j]
                            if (protLUT.has_key(substrate2)):
                                for peptide2 in protLUT[substrate2]:
                                    newEdge = network.NetworkEdge(peptide1, peptide2, score=1)
                                    edgeName = peptide1.strip() + '.' + peptide2.strip()
                                    if (not edgeList.has_key(edgeName)):
                                        edgeList[edgeName] = newEdge
            if (len(edgeList) > 0):
                networks.append(network.Network(name=factor, edges=edgeList, weight=score))
        return(networks)

class PTMFactor:
    """
    Store the known protein level PTM regulation 
    """
    def __init__(self, factor, substrate, coregulator="NA"):
        """
        Constructor
        """
        self._factor_ = factor
        self._substrate_ = substrate
        self._coregulator_ = coregulator
   
        
if __name__ == '__main__':
    print('Debugging ptmPeptide.py')
    test = PTMDataFactory()
    curDM = test.returnDM()
    curNets = test.returnNetworks()

    print('Writing testDM.tsv')
    curDM.write_tsv_file('./testDM.tsv')
    
