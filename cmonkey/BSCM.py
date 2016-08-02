"""BSCM.py - Module for Bicluster Sampled Coherence Matrix

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.

To Do: Write a 'save' and 'load' function that pickles/writes and unpickles/reads.
      Integrate that with main cMonkey save and load functions.

To Do: Write a function for resplitting clusters based on a new ratios matrix

This module implements the algorithm described in:
Bicluster Sampled Coherence Metric (BSCM) provides an accurate environmental context for phenotype predictions
Danziger et al.
BMC Systems Biology, 2015
"""
import numpy as np
import scipy as sp
import math
import random
import datetime as dt
import multiprocessing as mp
import logging

import cmonkey.datamatrix as datamatrix
import cmonkey.util as util


def getVarianceMeanSDvect_mp_wrapper(args):
    return getVarianceMeanSDvect(args[0], args[1], args[2], args[3], args[4], args[5], args[6])


def getVarianceMeanSDvect(ratioVect, n, tolerance = 0.01, maxTime=600, chunkSize=200, verbose=False,
                          expName=None):
    """Given a ratios matrix and a number of genes, figure out the expected distribution of variances
       Will sample background until the mean and sd converge or the operation times out
       Will return a list of variances to be used for statistical tests,
       or return nan if only nan values in ratioVect

     Keyword arguments:
     ratioVect  -- A a vector of ratios
     n          -- The number of genes to sample
     tolerance  -- The fraction tolance to use as a stopping condition (DEFAULT: 0.01)
     maxTime    -- The approximate maximum time to run in seconds (DEFAULT: 600)
     chunkSize  -- The number of samples to add between test (DEFAULT: 200)
     verbose    -- Set to false to suppress output (DEFAULT: False)
     expName    -- Set to echo this name if verbose = True (DEFAULT: None)

     Useage:
     varDist = getVarianceMeanSD(ratioVect, n)
    """
    ratioVect = [x for x in ratioVect if not math.isnan(x)]

    if verbose:
        logging.info("Calculating background for %d sampled from %d in %s", n, len(ratioVect), expName)

    if n <= 1 or n > len(ratioVect):
        return [np.nan]

    varList = []
    repeat = True
    startTime = dt.datetime.now()

    while repeat:
        newVars = []
        for i in range(0, chunkSize):
            curSample = random.sample(ratioVect, n)
            try:
                newVar = np.var(curSample)
            except:
                newVar = 0
            newVars.append(newVar)

        if len(varList) > 0: #True if past the first sample
            oldMean = np.mean(varList)
            oldVar = np.var(varList)
            varList = varList+newVars
            newMean = np.mean(varList)
            newVar = np.var(varList)
            meanWinTol = abs(newMean-oldMean) < tolerance*abs(oldMean)
            varWinTol = abs(oldVar-newVar) < tolerance*abs(oldVar)
            if meanWinTol and varWinTol:
                repeat = False
        else:
            varList = varList+newVars

        curTime = dt.datetime.now()
        if (curTime-startTime).seconds > maxTime:
            repeat = False

    return varList


class BSCM:
    """This is a class is designed to sample N items from a single vector
    until it reaches a certain convirgence criteria.  Once that's
    completed, it can be queried to return a p-Value for a specific set of genes
    Right now it copies ratios, which will waste some memory
    """
    def __init__(self, ratios, tolerance = 0.001, maxTime=600, chunkSize=200, verbose=False, useChi2=False):
        """Given a ratios matrix and a number of genes, figure out the expected distribution of variances
           Will sample background until the mean and sd converge or the operation times out

         Keyword arguments:
         ratios         -- A DataMatrix object from 'cmonkey.datamatrix'
         tolerance      -- The fraction tolance to use as a stopping condition (DEFAULT: 0.001)
         maxTime        -- The approximate maximum time to run in seconds (DEFAULT: 600)
         chunkSize      -- The number of samples to add between test (DEFAULT: 200)
         verbose        -- Set to false to suppress output (DEFAULT: False)
         useChi2        -- Set to True to fit a chi2 instead of storing all values (DEFAULT: False)
        """
        self.allVars = {} #Store all of the variances here.  Structure: allVars[expName][numExp]
        self.ratios = ratios
        self.tolerance = tolerance
        self.maxTime = maxTime
        self.chunkSize = chunkSize
        self.verbose = verbose
        self.useChi2 = useChi2

    def getPvals(self, geneNames, num_cores=1):
        """Get p-Values for the the list of genes, one for each column in the ratios matrix

         Keyword arguments:
         geneNames  -- A list of genes in the cluster
         singleCore -- Set to True to use a single core.  False may not work.
        """
        pVals = {}

        relGenes = list(set(geneNames) & set(self.ratios.row_names))
        curGeneMatrix = self.ratios.submatrix_by_rows(self.ratios.row_indexes_for(relGenes))

        noVarNs = []  #These three matrices should have a matched order
        noVarRats = []  #It would be better to have a single list
        noVarCns = []   #With 3 named elements in each list item

        for cn in self.ratios.column_names:
            colIdx = curGeneMatrix.column_indexes_for(column_names = [cn])
            geneVect = curGeneMatrix.column_values(column = colIdx)
            geneVect = [x for x in geneVect if not math.isnan(x)]
            n = len(geneVect)

            if not self.allVars.get(cn, False):
                self.allVars[cn] = {}

            #For loop: efficiently use multicore by precalculating additional numbers of genes
            i_s = [n]
            if num_cores > 1:
                i_s = [n-3, n-2, n-1, n, n+1, n+2, n+3]
            for i in i_s:
                if not self.allVars[cn].get(str(i), False) and i >= 0:
                    ratioVect = self.ratios.column_values(column = self.ratios.column_indexes_for(column_names = [cn]))
                    noVarNs.append(i)
                    noVarRats.append(ratioVect.tolist())
                    noVarCns.append(cn)

        #  2) Use a pool of workers to calculate a distribution for each of the tuples
        if len(noVarNs) > 0:
            logging.info("Calculating some backgrounds for about %d genes", len(geneNames))
            if self.useChi2:
                logging.info("\tFitting variance samples to Chi2 distribution")

            if num_cores > 1:
                newargs = []
                for i in range(0, len(noVarNs)):
                    newargs.append([noVarRats[i], noVarNs[i], self.tolerance,
                                    self.maxTime, self.chunkSize, self.verbose, noVarCns[i]])
                pool = mp.Pool(num_cores)
                newVars = pool.map(getVarianceMeanSDvect_mp_wrapper, newargs)
                pool.close()
                pool.join()
            else:
                tolerance = np.repeat(self.tolerance, len(noVarNs)).tolist()
                maxTime = np.repeat(self.maxTime, len(noVarNs)).tolist()
                chunkSize = np.repeat(self.chunkSize, len(noVarNs)).tolist()
                verbose = np.repeat(self.verbose, len(noVarNs)).tolist()
                newVars = list(map(getVarianceMeanSDvect, noVarRats, noVarNs, tolerance, maxTime,
                                   chunkSize, verbose, noVarCns))

        #  3) Assign the new values into the empty slots
        for idx in range(0,len(noVarCns)):
            cn = noVarCns[idx]
            curN = str(noVarNs[idx])
            if self.useChi2:
                curVars = newVars[idx]
                self.allVars[cn][curN] = sp.stats.chi2.fit(curVars, df=int(curN))
            else:
                self.allVars[cn][curN] = newVars[idx]

        #  4) Calculate the p-Values
        pVals = {}
        for cn in self.ratios.column_names:
            colIdx = curGeneMatrix.column_indexes_for(column_names = [cn])
            geneVect = curGeneMatrix.column_values(column = colIdx)
            geneVect = [x for x in geneVect if not math.isnan(x)]
            n = str(len(geneVect))

            if len(geneVect) <= 1 or np.any(np.isnan(self.allVars[cn][str(n)])):
                pVals[cn] = 1
            else:
                curVar = np.var(geneVect)
                if self.useChi2:
                    [df, loc, scale] = self.allVars[cn][str(n)]
                    pVals[cn] = 1-sp.stats.chi2.sf(curVar, df=df, loc=loc, scale=scale)
                else:
                    pVals[cn] = np.mean(self.allVars[cn][str(n)] < curVar)

        return pVals

    def resplit_clusters(self, membership, cutoff=0.05):
        """Get p-Values for the the list of genes, one for each column in the ratios matrix
        Note: this will increase the number of elements in each row of 'membership.col_membs'

         Keyword arguments:
         membership -- A membership object containing all of the cluster membership information
         cutoff     -- The p-Value inclusion cutoff (DEFAULT: 0.05)
        """
        pDict = {}
        for cluster in range(1, membership.num_clusters() + 1):
            cur_genes = membership.rows_for_cluster(cluster)
            cur_pvals = self.getPvals(geneNames=cur_genes, num_cores=1)

            for curCol in cur_pvals.keys():
                if not curCol in pDict:
                    pDict[curCol] = []

                if cur_pvals[curCol] <= cutoff:
                    pDict[curCol].append(cluster)
                else:
                    pDict[curCol].append(0)

        membership.col_membs = np.zeros((len(membership.col_membs),membership.num_clusters()),
                                        dtype='int32')
        for col in pDict.keys():
            membership.col_membs[membership.colidx[col]] = np.array(pDict[col], dtype='int32')

        return membership
