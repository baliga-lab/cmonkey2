#!/usr/bin/env python
# vi: sw=4 ts=4 et:
"""make.sequences.py
Creates sequences file from genome fasta and gene model
feature file. Use as script. Only needs to be run manually once per config,
instead of each time cMonkey runs. Then supply sequences file to
GenericOrganism.

Default is to count coordinates from the gene model start. Can also output 3' 'UTR' sequences using the '--fromend' option.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""

import string
import seqtools as st
from optparse import OptionParser
import util

def make_sequences( genome_fasta_file, gene_features_file,
        outfile='sequences.csv', distance={'upstream':300,'downstream':100}, from_end=False, fasta=False ):

    if from_end:
        distance = ( distance['upstream'], distance['downstream'] )
    else:
        '''WARNING: as of 2012-03-22, the st.extract functions used flipped distances!
           e.g. distance[1] is the UPSTREAM distance and distance[0] is the DOWNSTREAM
           CHECK YOUR SEQUENCES after running this! Also, a negative number is expected for
           DOWNSTREAM. So, (-100,300) must be passed to st.extract_upstream in order to get
           a sequence from 300 upstream to 100 downstream. WEIRD!'''
        distance = (-1*distance['downstream'],distance['upstream'])

    contig_sequences = st.read_sequences_from_fasta_file( genome_fasta_file )
    # convert contig_sequences to dictionary (this func returns a list of tuples)
    contig_dict = {}
    for name, seq in contig_sequences:
        contig_dict[name] = seq
    print 'loaded %i contigs' %len(contig_dict)
    print string.join( [ '%s: %ibp' %(a,len(b)) for a,b in contig_dict.items()] , ',' )

    features = st.read_features_from_file( gene_features_file )
    print 'loaded %i features' %len(features)
#    print str(features.values()[1])

    sequences = []
    for feature in features.values():
        location = feature.location
#        print location, location.contig, distance, feature.id()
        if from_end:
            sequences.append( ( feature.id(), st.extract_downstream(contig_dict[location.contig], location, distance)[1] ) )
        else:
            sequences.append( ( feature.id(), st.extract_upstream(contig_dict[location.contig], location, distance)[1] ) )
#        print sequences[feature.id()]

    outf = open(outfile,'w')
    if fasta: st.write_sequences_to_fasta_file(outf,sequences)
    else:
        sep = ','
        for id, seq in sequences:
            outf.write( '%s%s%s\n' %(id,sep,seq) )
    outf.close()

if __name__ == "__main__":
    p = OptionParser()
    p.add_option('-o','--output',default='sequences.csv')
    p.add_option('-u','--upstream',type='int',default=250)
    p.add_option('-d','--downstream',type='int',default=50)
    p.add_option('-e','--fromend',action='store_true',default=False,help='relative to gene end instead of gene start')
    p.add_option('-f','--fasta',action='store_true',default=False,help='output in fasta format')
    opt,args = p.parse_args()

    if len(args) < 2:
        print('Usage: ./make.sequences.py <genome-fasta-file> <gene-features-file> [options]')

    make_sequences( args[0], args[1], opt.output, distance={'upstream':opt.upstream, 'downstream':opt.downstream}, from_end=opt.fromend, opt.fasta )
    print 'NOW VERIFY YOUR SEQUENCES! There are differing formats for supplying upstream and downstream shifts! (e.g. (200,100), (-200,100), (100,200)....) This script assumes that: -u 200 -d 100 means: from 200 upstream to 100 downstream'
