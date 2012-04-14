'''
Created on Mar 25, 2012

@author: Frank Schmitz, SBRI, 2012 -
'''


MIN_SEQ_LENGTH = 8      # minimum Sequence Length to be included in search

import cPickle
import logging
import organism as org
import thesaurus as th
import util
from Einheit import *

class Kapsel():
    '''
    a class that contains a necc stuff to generate a
    a microbe or whatever organism, preferrably an eucaryote
    or in particular mouse or human
    from a directory with RSAT files
    will generate a thesaurus and so on...
    '''


    def __init__(self, org_code, Kapsel_name, Kapsel_description):
        '''
        Constructor with the very bare bones
        everything else is implemented within the specific methods
        '''
        self.Kapsel_organism = org_code
        self.Kapsel_name = Kapsel_name
        self.Kapsel_description = Kapsel_description


        # something borrowed...
        self.kegg_organism = ""
        self.__rsat_info = ""
        self.__microbes_online_db = ""
        self.go_taxonomy_id = ""
        self.__search_distances = ""
        self.__scan_distances = ""

        self.__synonyms = {}
        self.TRex = None
        self.__networks = None          #default for empty
        self.__operon_mappings = None  # lazy loaded
        self.__network_factories = []

        self.__GeneDict = {}
        self.__contigDict = {}

        self.__sequenceMethod = None 
        
        logging.info("\x1b[31mKapsel:\t\x1b[0mcreated: %s" %(self.Kapsel_description))

    
    def create(self):
        pass
    
    def setSequenceMethod(self, method, filetype):

        self.__sequenceMethod = method
        self.__sequenceFileType = filetype

    def getSequenceMethod(self):
        if self.__sequenceMethod == None:
            # override with default
            self.__sequenceMethod = "de novo"   # get it from the huuuuuge files.. 
        return self.__sequenceMethod

    def getSequenceFileType(self):
        if self.__sequenceFileType == None:
            # override with default
            self.__sequenceFileType = "repeat_masked"   # get it from the repeat Masked File 
        return self.__sequenceFileType
        
    
    def addNWfactory(self, factory, weight, nwfactname = "no Name"):
        # this one adds a NW factory to the list, to ve used later on
        self.__network_factories.append((factory, weight))
        logging.info("\x1b[31mKapsel:\t\x1b[0mNW factory %s added" %(nwfactname))

    def networks(self):
        """returns this organism's networks"""
        if self.__networks == None:
            self.__networks = []
            for make_network, weight in self.__network_factories:
                self.__networks.append(make_network(self))
            for NW in self.__networks:
                NW.TranslateNetwork(self.thesaurus())
        return self.__networks

    def getContig(self):
        pass

    
    def addThesaurusFile(self, thFile, thFileType):
        # here we add a thesaurus file, can choose from different filetypes
        # at least later
        if thFileType == "RSAT":
            infile = util.DelimitedFile.read(thFile, sep='\t', has_header=False, comment = "--")
            tempD = th.create_from_rsat_feature_names(infile)
            self.__synonyms = self.MergeDicts(self.__synonyms, tempD)
        else:
            thFileType == "unassigned - ERROR"
        logging.info("\x1b[31mKapsel:\t\x1b[0mThesaurus added, %s type" % (thFileType))
        
    def thesaurus(self):
        """reads the thesaurus from a feature_names file. The thesaurus
        is also cached, because it is used many times
        """
        if not self.__synonyms:
            logging.warning("no thesaurus loaded!")
        return self.__synonyms
    def put_thesaurus(self, Dict):
        """
        dominates the thesaurus - if given at all
        """
        self.__synonyms = Dict
    def TRex(self):
        """reads the thesaurus from a feature_names file. The thesaurus
        is also cached, because it is used many times
        """
        if not self.__TRex:
            logging.warning("Kapsel:\tno TRex loaded loaded!")
        return self.__TRex
    def put_TRex(self, Dict):
        """
        dominates the thesaurus - if given at all
        """
        self.__TRex = Dict

    def put_ReSaurus(self, Dict):
        """
        define the ReSaurus - the Dict re-translating for reporting and plotting
        
        """
        self.__ReSaurus = Dict
    def get_ReSaurus(self):
        return self.__ReSaurus


    def MergeDicts(self, Dictone, Dicttwo):
        # a lame function to merge two dictionaries
        # allows for multiple thesaurus files
        
        for k, v in Dicttwo.iteritems():
            if k not in Dictone:
                Dictone[k] = v
        return Dictone
        
    def addRSATenvir(self, featdir, featfile, UTRfile):
        """
        adds the environmental information for the RSAT dir and the features
        only local loading ins implemented, only for eucaryotes
        """
        self.__RSATDir = featdir
        self.__featfile = featfile
        self.__UTRfile = UTRfile
    
    
    def addRSAT_feature(self, grabtype, Sboundaries, Bboundaries):
        """
        this method slurps the features / cds file from RSAT into a collection of
        gene - specific sequence / subunits
         
        """
        boundaries = Sboundaries + Bboundaries
        infile = util.DelimitedFile.read(self.__featfile, sep='\t', has_header=False, comment = "--")
        logging.info("\x1b[31mKapsel:\t\x1b[0mread features/cds File, type %s" % (grabtype))
        for line in infile.lines():
            geneID = line[9]
            geneTransID = line[9]
            geneProduct = line[0]
            description = line[7]
            contig = line[3]
            Start = int(line[4])
            Stop = int(line[5])
            strand = line[6]
            transcript = line[8]
            # generate Gene instance
            try:
                new_geneID = self.thesaurus()[geneTransID]
            except KeyError:
                try:
                    new_geneID = self.thesaurus()[geneProduct]
                except KeyError:
                    try:
                        new_geneID = self.thesaurus()[geneID]
                    except KeyError:
                        new_geneID = geneID
            try:
                gene = self.__GeneDict[new_geneID]
            except KeyError:
                gene = Einheit(geneProduct, description, new_geneID)
                # set contig
                gene.setContig(contig)
            subUnit = SubEinheit(transcript)
            subUnit.pushSeq(contig, Start, Stop, strand, "Promoter", boundaries)
            gene.pushSubUnit(subUnit)
            # put into list of contig contents
            try:
                tempL = self.__contigDict[contig]
            except KeyError:
                tempL = []
            if new_geneID not in tempL:
                tempL.append(new_geneID)
            self.__contigDict[contig] = tempL

            # save gene into Dict    
            self.__GeneDict[new_geneID] = gene
    
    def addRSAT_UTR(self, grabtype, Sboundaries, Bboundaries):
        """
        this method slurps the features / cds file from RSAT into a collection of
        gene - specific sequence / subunits
         

        """
        boundaries = Sboundaries + Bboundaries
        infile = util.DelimitedFile.read(self.__UTRfile, sep='\t', has_header=False, comment = "--")
        logging.info("\x1b[31mKapsel:\t\x1b[0mread UTR File, type %s" % grabtype)
        for line in infile.lines():
            if line[1] != grabtype:
                continue
            geneProduct = line[0]   # well, its just a UTR
            geneID = line[9]
            transcript = line[8]
            contig = line[3]
            Start = int(line[4])
            Stop = int(line[5])
            strand = line[6]
            description = line[7]
            # generate Gene instance - only if needed!!!
 
            try:
                new_geneID = self.thesaurus()[transcript]
            except KeyError:
                try:
                    new_geneID = self.thesaurus()[geneID]
                except KeyError:
                    new_geneID = geneID

            try:
                gene = self.__GeneDict[new_geneID]
            except KeyError:
                gene = Einheit(geneProduct, description, new_geneID)
                # set contig
                gene.setContig(contig)

            subUnit = SubEinheit(transcript)
            subUnit.pushSeq(contig, Start, Stop, strand, "3pUTR", boundaries)
            
            gene.pushSubUnit(subUnit)


            # put into list of contig contents
            try:
                tempL = self.__contigDict[contig]
            except KeyError:
                tempL = []
            if new_geneID not in tempL:
                tempL.append(new_geneID)
            self.__contigDict[contig] = tempL

            # save gene into Dict    
            self.__GeneDict[new_geneID] = gene

    
    def CollectSequences(self, type):
        logging.info("\x1b[31mKapsel:\t\x1b[0mcollecting Sequences")
        contigs = sorted(self.__contigDict.keys())
        for contig in contigs:
            content = self.getContigContent(contig, type)
            logging.info("\x1b[31mKapsel:\t\x1b[0mparsing contig file %s, %s type" % (contig, type))
            for geneID in self.__contigDict[contig]:
                gene = self.__GeneDict[geneID]
                gene.updateEinheit_Seq(content)




    def getContigContent(self, contig, contigType):
        contigname = contig.replace(":", "_")
        if contigType == "repeat_masked":
            contigname = contigname + "_repeat_masked"
        contigname = contigname + ".raw"
        bigfile = self.__RSATDir + contigname
        handle = open(bigfile, "r")
        content = handle.read()
        handle.close
        return content
    

    def getUnit(self, GeneName):
        try:
            return self.__GeneDict[GeneName]
        except KeyError:
            logging.debug("gene not found")
            return Einheit("NullGene", "DefaultReturn for non found Gene", "ID0000")
            

    def sequences_for_genes_scan(self, gene_aliases, seqtype):
        newAliases = []
        for gene in gene_aliases:
            try:
                gA = self.__synonyms[gene]
                newAliases.append(gA)
            except KeyError:
                #print "%s not found in thesaurus - fuck!" % (gene)
                pass
        logging.info("working on %s genes for background model generation" % (len(newAliases)))
        tempD = {}
        if seqtype == None:
            return {}
        for Unit in self.__GeneDict.values():
            arbicount = 0
            if Unit.getGeneID() not in newAliases and Unit.getGeneName() not in newAliases: continue
            BSeqD = Unit.getSeqDs()[1]
            
            try:
                SeqL = BSeqD[seqtype]
            except KeyError:
                SeqL = []
            for Seq in SeqL:
                tempD['%s_%03d' % (Unit.getGeneID(), arbicount)] = ("-", Seq)
                arbicount += 1
        logging.info("collected %s sequences for background model generation" % (len(tempD)))
        return tempD


    def sequences_for_genes_search(self, gene_aliases, seqtype):
        newAliases = []
        for gene in gene_aliases:
            try:
                gA = self.__synonyms[gene]
                newAliases.append(gA)
            except KeyError:
                pass
        #logging.info("working on %s genes for model generation" % (len(newAliases)))
        tempD = {}
        if seqtype == None:
            return {}
        for Unit in self.__GeneDict.values():
            arbicount = 0
            if Unit.getGeneID() not in newAliases and Unit.getGeneName() not in newAliases: continue
            SeqD = Unit.getSeqDs()[0]
            
            try:
                SeqL = SeqD[seqtype]
            except KeyError:
                SeqL = []
            for Seq in SeqL:
                if len(Seq) < MIN_SEQ_LENGTH: continue
                tempD['%s_%03d' % (Unit.getGeneID(), arbicount)] = Seq
                arbicount += 1
        #logging.info("collected %s sequences for motif generation" % (len(tempD)))
        return tempD

        
    def feature_ids_for(self, gene_aliases):
        """Helper method to retrieve a list of feature_ids for the
        specified alias list"""
        synonyms = self.thesaurus()
        return [synonyms[alias] for alias in gene_aliases if alias in synonyms]
   
    def is_eukaryote(self):
        """Determines whether this object is an eukaryote"""
        return True


    def __sequences_for_genes(self, seqtype, gene_aliases, distance):
        """retrieves the specified sequences from the supplied genomic data"""
        if not seqtype in self.__seqs:
            dfile = util.DelimitedFile.read(self.__seq_filenames[seqtype], sep=',')
            self.__seqs[seqtype] = {}
            for line in dfile.lines():
                self.__seqs[seqtype][line[0].upper()] = line[1]
        result = {}
        for alias in gene_aliases:
            if alias in self.thesaurus():
                gene = self.thesaurus()[alias]
                if gene in self.__seqs[seqtype]:
                    result[gene] = self.__seqs[seqtype][gene]
                else:
                    #logging.warn("Gene '%s' not found in 3' UTRs", gene)
                    pass
            else:
                #logging.warn("Alias '%s' not in thesaurus !", alias)
                pass
        return result

        
    
