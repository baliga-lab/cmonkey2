#!/usr/bin/env python3
"""all_tests.py - run all unit tests in the project

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import unittest
import xmlrunner

import orig_membership_test as omembtest
import datamatrix_test as dmtest
import util_test as ut
import organism_test as ot
import seqtools_test as stt
import thesaurus_test as tht
import operon_nw_test as opnwt
import network_test as nwt
import microarray_test as mat
import meme_test as met
import pssm_test as pt
import combiner_test as ct
import read_wee_test as rwt

import rsat_test
import meme430_test
import halo_genes_test
import microbes_online_test as mo_test
import iteration_test
import postproc_test
import setenrichment_test as se_test

# pylint: disable-msg=C0301
if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(dmtest.DataMatrixTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(dmtest.DataMatrixFactoryTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(dmtest.NoChangeFilterTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(dmtest.CenterScaleFilterTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(dmtest.QuantileNormalizeTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ut.DelimitedFileTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ut.UtilsTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ut.LevenshteinDistanceTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ut.BestMatchingLinksTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ut.Order2StringTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ot.MicrobeTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(stt.SeqtoolsTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(stt.FastaTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(stt.LocationTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(tht.DelimitedFileFactoryTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(opnwt.ReadOperonNetworkTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(opnwt.GetOperonPairsTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(nwt.NetworkTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(omembtest.OrigMembershipTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(mat.ComputeArrayScoresTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(met.MemeTest))

    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(pt.PssmTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(ct.CombinerTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(rwt.ReadWeeTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(se_test.DiscreteEnrichmentSetTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(se_test.CutoffEnrichmentSetTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(se_test.SetTypeTest))

    # web based tests
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
        rsat_test.RsatDatabaseTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
        mo_test.MicrobesOnlineTest))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
        meme430_test.Meme430Test))
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
        halo_genes_test.HaloGeneTest))

    # iteration_test
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
        iteration_test.IterationTest))

    # postproc_test
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(
        postproc_test.PostprocTest))

    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
      xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(SUITE))
    else:
      unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
