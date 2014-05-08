"""extract_string_links.py - Helper tool for generating a preprocessed
STRING link file for a specific organism

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import re
import util
import stringdb
import rsat
import organism
import thesaurus
import gzip
import network


# This tool extracts the relevant protein-protein interaction links for a specific
# organism from a huge STRING database file. Since this file contains protein links
# for a lot of different organisms, processing has to be performed to be practicable.
# The gene links are taken that are found in the RSAT database for the specified
# organism
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
CACHE_DIR = 'cache'
KEGG_FILE = 'testdata/KEGG_taxonomy'


def process_stringdb(stringfile, taxonomy_id):
    result = []
    max_score = 0.0
    line = stringfile.readline()  # ignore header
    print "processing STRING database..."
    line = stringfile.readline()
    prefix = "%s." % str(taxonomy_id)
    while line != None:
        if line.startswith(prefix):
            print line.strip()
        line = stringfile.readline()

"""
    while line != None:
        comps = line.strip().split(' ')
        if len(comps) >= 3:
            gene1 = re.sub('^\d+[.]', '', comps[0])
            gene2 = re.sub('^\d+[.]', '', comps[1])
            score = float(comps[2])
            if gene1 in synonyms and gene2 in synonyms:
                if score > max_score:
                    max_score = score
                result.append(network.NetworkEdge(gene1, gene2, score))
        line = stringfile.readline()

    print "processed STRING database, normalizing..."
    normalized_result = stringdb.normalize_edge_list(result, max_score)
    print "edges normalized."
    for edge in normalized_result:
        print "%s\t%s\%f" % (edge.source(), edge.target(), edge.score())
"""

if __name__ == '__main__':
    print("extract_string_links.py, (c) 2012, Institute for Systems Biology")
    print('This program is licensed under the General Public License V3.')
    print('See README and LICENSE for details.\n')
    if len(sys.argv) <= 2:
        print('Usage: python extract_string_links.py <stringdb-path> <organism-code>')
    else:
        rsatdb = rsat.RsatDatabase(RSAT_BASE_URL, CACHE_DIR)
        kegg_mapper = organism.make_kegg_code_mapper(util.DelimitedFile.read(KEGG_FILE, sep='\t',
                                                                             has_header=True,
                                                                             comment='#'))
        kegg_org = kegg_mapper(sys.argv[2])
        rsat_info = organism.RsatSpeciesInfo(rsatdb, kegg_org, None, None)
        print "RSAT SPECIES: ", rsat_info.species
        print "TAX ID: ", rsat_info.taxonomy_id
        feature_names = rsatdb.get_feature_names(rsat_info.species)
        feature_names_dfile = util.DelimitedFile.create_from_text(feature_names, comment='--')
        synonyms = thesaurus.create_from_rsat_feature_names(feature_names_dfile)
        string_filepath = sys.argv[1]
        if string_filepath.endswith('.gz'):
            with gzip.open(string_filepath) as stringfile:
                process_stringdb(stringfile, rsat_info.taxonomy_id)
        else:
            with open(string_filepath) as stringfile:
                process_stringdb(stringfile, rsat_info.taxonomy_id)
