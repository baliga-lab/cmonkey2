#!/usr/bin/python
import argparse
import psycopg2
import os
import os.path
import util
import organism as org
import addnwportal

PSQL = "dbname='network_portal' user='dj_ango' host='localhost' password='django'"
BASEURL = 'http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId='  #224308;export=tab

def mark_regulators(pgconn, species_id, regulators, synonyms):
    print "marking %d transcription factors..." % len(regulators)
    cursor = pgconn.cursor()
    for regulator in regulators:
        primary = synonyms.get(regulator, regulator)
        cursor.execute('update networks_gene set transcription_factor = true where species_id = %s and name = %s', (species_id, primary))
    pass

if __name__ == '__main__':
    description = "fix_genes.py - fixup network portal genes"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--organism', required=True, help='KEGG organism code')
    args = parser.parse_args()

    tffile = "%s_TF.txt" % args.organism
    with open(os.path.join('pipeline_organisms/regulators', tffile)) as infile:
        regulators = [line.strip() for line in infile]

    pgconn = psycopg2.connect(PSQL)
    pgconn.set_isolation_level(0)  # set auto commit
    cursor = pgconn.cursor()
    cursor.execute('select id, ncbi_taxonomy_id from networks_species where short_name = %s',
                   [args.organism])
    species_id, taxonomy_id = cursor.fetchone()
    cache_filename = 'cache/ncbi_' + str(taxonomy_id)
    url = BASEURL + ("%d;export=tab" % taxonomy_id)

    microbedb, organism = addnwportal.make_microbe(args.organism)
    synonyms = organism.thesaurus()

    mark_regulators(pgconn, species_id, regulators, synonyms)  # mark transcription factors
    util.get_url_cached(url, cache_filename)
    with open(cache_filename) as cached_file:
        lines = [line.split('\t') for i, line in enumerate(cached_file) if i > 0]

    found = 0
    not_found = 0
    for i, line in enumerate(lines):
        sysname = line[7]  # mapped to name
        name = line[8]     # mapped to common name
        description = line[9]
        #accession = line[1]
        primary = None
        tomap = None
        if sysname in synonyms:
            primary = synonyms[sysname]
            tomap = sysname
        elif name in synonyms:
            primary = synonyms[name]
            tomap = name

        if primary != None:
            #print "%s -> %s" % (tomap, primary)
            cursor.execute('select id from networks_gene where name = %s', [primary])
            result = cursor.fetchall()
            if len(result) > 0:
                cursor.execute('update networks_gene set name = %s, common_name = %s, description = %s where id = %s', [sysname, name, description, result[0][0]])
                #print "remapping '%s'" % primary
                found += 1
            else:
                #print "not found '%s'" % primary
                not_found += 1
    print "remapped %d of %d possible." % (found, found + not_found)
