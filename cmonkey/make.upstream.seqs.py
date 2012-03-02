#!/usr/bin/env python
# vi: sw=4 ts=4 et:
"""make.upstream.sequences.py
Creates upstream sequences file from genome fasta and gene model
feature file. Use as script. Only needs to be run manually once per config,
instead of each time cMonkey runs. Then supply sequences file to
GenericOrganism.

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""

import string
import seqtools as st
from optparse import OptionParser
import util

def make_upstream_sequences( genome_fasta_file, gene_features_file,
        outfile='upstream.sequences', distance=[-300,100] ):
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
        location = feature.location()
#        print location, location.contig, distance, feature.id()
        sequences.append( ( feature.id(), st.extract_upstream(contig_dict[location.contig], location, distance)[1] ) )
#        print sequences[feature.id()]

    outf = open(outfile,'w')
#    st.write_sequences_to_fasta_file(outf,sequences)
    sep = ','
    for id, seq in sequences:
        outf.write( '%s%s%s\n' %(id,sep,seq) )
    outf.close()

if __name__ == "__main__":
    p = OptionParser()
    p.add_option('-o','--output',default='upstream.sequences')
    p.add_option('-u','--upstream',type='int',default=-250)
    p.add_option('-d','--downstream',type='int',default=50)
    opt,args = p.parse_args()

    if len(args) < 2:
        print('Usage: ./make.upstream.sequences.py <genome-fasta-file> <gene-features-file> [options]')

    make_upstream_sequences( args[0], args[1], opt.output, distance=[opt.upstream, opt.downstream] )
