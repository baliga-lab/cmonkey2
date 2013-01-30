#!/usr/bin/python
import argparse
import sqlite3
import organism as org
import util, rsat, microbes_online
import psycopg2
import collections
import datamatrix as dm

# Script to add a cMonkey python run to the network portal database
KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
STRING_URL_PATTERN = "http://networks.systemsbiology.net/string9/%s.gz"
CACHE_DIR = 'cache'
# Map for UCSC genome browser
UCSC_MAP = {'gsu': 'geobSulf', 'cac': '', 'bce': '', 'bsu': 'baciSubt2',
            'rsp': 'rhodSpha', 'cje': 'campJeju', 'ype': 'yersPest_CO92',
            'bth': 'bactThet_VPI_5482', 'ttj': 'therTher_HB8',
            'pae': 'pseuAeru'
             }

MicrobeDB = collections.namedtuple('MicrobeDB', ['keggfile', 'rsatdb', 'rsat_info'])

def make_microbe(code):
    """assemble organism related information and return it to the caller"""
    keggfile = util.read_dfile(KEGG_FILE_PATH, comment='#')
    rsatdb = rsat.RsatDatabase(RSAT_BASE_URL, CACHE_DIR)
    kegg_mapper = org.make_kegg_code_mapper(keggfile)
    rsat_mapper = org.make_rsat_organism_mapper(rsatdb)
    rsat_info = rsat_mapper(kegg_mapper(code))
    microbedb = MicrobeDB(keggfile, rsatdb, rsat_info)
    print "NCBI CODE IS: ", rsat_info.taxonomy_id
    gofile = util.read_dfile(GO_FILE_PATH)
    mo_db = microbes_online.MicrobesOnline()
    search_distances= {'upstream': (-20, 150)}
    scan_distances = {'upstream': (-30, 250)}
    org_factory = org.MicrobeFactory(kegg_mapper,
                                     rsat_mapper,
                                     org.make_go_taxonomy_mapper(gofile),
                                     mo_db,
                                     [])
    organism = org_factory.create(code, search_distances, scan_distances)
    return microbedb, organism

def add_species(pgconn, orgcode, species, ncbi_code, ucsc_code):
    """add a species entry for the specified organism"""
    cur = pgconn.cursor()
    cur.execute('select id from networks_species where short_name = %s', [orgcode])
    result = cur.fetchall()
    if len(result) == 0:
        print "adding species: ", orgcode
        cur.execute('insert into networks_species (name, short_name, ncbi_taxonomy_id, ucsc_id, created_at) values (%s, %s, %s, %s, now()) returning id', [species, orgcode, ncbi_code, ucsc_code])
        species_id = cur.fetchone()[0]
    else:
        print "retrieving species: ", orgcode
        species_id = result[0][0]

    cur.close()
    return species_id

def add_network(pgconn, species_id, species):
    """add a network entry for the specified organism"""
    cur = pgconn.cursor()
    cur.execute('select id from networks_network where species_id = %s', [species_id])
    result = cur.fetchall()
    if len(result) > 0:
        print "retrieving network for species: ", species
        return result[0][0]
    else:
        print "creating network for species: ", species
        name = '%s network' % species
        datasource = 'MicrobesOnline & cMonkey/Python'
        description = 'regulatory network of %s from cMonkey/Python pipeline' % species
        print "inserting new network '%s'" % name
        cur.execute("insert into networks_network (species_id, name, data_source, description, created_at) values (%s, %s, %s, %s, now()) returning id", [species_id, name, datasource, description])
        return cur.fetchone()[0]

def add_biclusters(pgconn, sqliteconn, network_id, iteration):
    """add or retrieve the biclusters specified in the cmonkey run"""
    cur = pgconn.cursor()
    cur.execute('select id, k, residual from networks_bicluster where network_id = %s',
                [network_id])
    result = cur.fetchall()
    if len(result) > 0:
        print "retrieving biclusters..."
        return [(row[0], row[1], row[2]) for row in result]
    else:
        print "adding biclusters..."
        sqlite_cur = sqliteconn.cursor()
        sqlite_cur.execute('select cluster, residual from cluster_residuals where iteration = ?', [iteration])
        resids = [(row[0], row[1]) for row in sqlite_cur.fetchall()]
        result = []
        for cluster, residual in resids:
            cur.execute('insert into networks_bicluster (network_id, k, residual) values (%s, %s, %s) returning id',
                        [network_id, cluster, residual])
            result.append((cur.fetchone()[0], cluster, residual))
        return result

def add_chromosome(pgconn, species_id, contig):
    """add the specified chromosome to the database or return its entry"""
    cur = pgconn.cursor()
    cur.execute('select id from networks_chromosome where refseq = %s', [contig])
    result = cur.fetchall()
    if len(result) == 0:
        print 'adding chromosome'
        cur.execute("insert into networks_chromosome (species_id, name, length, topology, refseq) values (%s, %s, 0, 'circular', %s) returning id", [species_id, contig, contig])
        return cur.fetchone()[0]
    else:
        print 'retrieving chromosome'
        return result[0][0]

def add_gene(pgconn, species_id, chr_id, name, common_name, gene_type, start, end, strand):
    """add the specified gene to the database or return it's entry in the
    network portal database"""
    try:
        start = int(start)
        end = int(end)
    except:
        # in case start and end are not defined
        start = 0
        end = 0
    strand = '-' if strand == 'R' else '+'
    cur = pgconn.cursor()
    cur.execute("select id from networks_gene where species_id = %s and chromosome_id = %s and name = %s", [species_id, chr_id, name])
    result = cur.fetchall()
    if len(result) == 0:
        print "adding %s %s %s %d %d %s" % (name, common_name, gene_type, start, end, strand)
        cur.execute('insert into networks_gene (species_id, chromosome_id, name, common_name, type, start, "end", strand) values (%s, %s, %s, %s, %s, %s, %s, %s) returning id',
                    [species_id, chr_id, name, common_name, gene_type, start, end, strand])
        return (cur.fetchone()[0], name)
    else:
        return (result[0][0], name)

def add_rsat_genes(pgconn, species_id, microbedb):
    """Add the organism's genes that are found in the features file in RSAT"""
    print "adding genes from RSAT..."
    rows = microbedb.rsatdb.get_features(microbedb.rsat_info.species).split('\n')
    rows = [row.split('\t') for row in rows
            if not row.startswith('--') and len(row.strip()) > 0]  # ignore comments
    rows = [row for row in rows if row[10].isdigit()]  # only allow valid genes
    contigs = set([row[3] for row in rows])
    chromosomes = [(add_chromosome(pgconn, species_id, contig), contig) for contig in contigs]
    chr_map = {chrom[1]: chrom[0] for chrom in chromosomes}
    genes = []
    for row in rows:
        genes.append(add_gene(pgconn, species_id, chr_map[row[3]],
                              row[0], row[2], row[1], row[4], row[5], row[6]))
    return chr_map, genes


def add_ratio_genes(pgconn, species_id, ratios, thesaurus, chr_map, gene_map):
    """add genes to the database that are not found in RSAT, we need to reference
    them in some way"""
    print "checking ratio matrix for missing genes..."
    non_mapped = [gene for gene in ratios.row_names if gene not in thesaurus]
    chr_id = chr_map.items()[0][1]
    return [add_gene(pgconn, species_id, chr_id, gene, gene, 'DUMMY', 0, 0, '+')
            for gene in non_mapped]


def add_condition(pgconn, network_id, condition):
    """adds or retrieves the specified condition"""
    cur = pgconn.cursor()
    cur.execute('select id from networks_condition where network_id = %s and name = %s',
                [network_id, condition])
    result = cur.fetchall()
    if len(result) == 0:
        print "adding condition '%s'" % condition
        cur.execute('insert into networks_condition (network_id, name) values (%s, %s) returning id', [network_id, condition])
        return (cur.fetchone()[0], condition)
    else:        
        return (result[0][0], condition)

def add_conditions(pgconn, network_id, ratios):
    """add the conditions contained in the ratios matrix"""
    print "adding conditions from ratios matrix"""
    conditions = ratios.column_names
    return [add_condition(pgconn, network_id, cond) for cond in ratios.column_names]

def add_bicluster_gene(pgconn, cluster_id, gene_id):
    cur = pgconn.cursor()
    cur.execute('select * from networks_bicluster_genes where bicluster_id = %s and gene_id = %s', [cluster_id, gene_id])
    result = cur.fetchall()
    if len(result) == 0:
        cur.execute('insert into networks_bicluster_genes (bicluster_id, gene_id) values (%s, %s)',
                    [cluster_id, gene_id])
        return True
    else:
        return False

def add_bicluster_genes(pgconn, sqliteconn, thesaurus, bicl_map, gene_map, iteration):
    """copy biclusters from sqlite to postgres if not there"""
    print "copying bicluster-gene memberships..."
    clusters = sorted(bicl_map.keys())
    sqlcur = sqliteconn.cursor()
    for cluster in clusters:
        sqlcur.execute('select name from row_members m join row_names n on m.order_num = n.order_num where iteration = ? and cluster = ?', [iteration, cluster])
        cluster_id = bicl_map[cluster]
        genes = [row[0] for row in sqlcur.fetchall()]
        mgenes = []
        # use the thesaurus to get to the normalized genes
        for gene in genes:
            if gene in thesaurus:
                mgenes.append(thesaurus[gene])
            else:
                mgenes.append(gene)
        gene_ids = [gene_map[gene] for gene in mgenes]
        added = 0
        for gene_id in gene_ids:
            added += 1 if add_bicluster_gene(pgconn, cluster_id, gene_id) else 0
        print "genes added to cluster ", cluster, ": ", added


def add_bicluster_condition(pgconn, cluster_id, cond_id):
    cur = pgconn.cursor()
    cur.execute('select * from networks_bicluster_conditions where bicluster_id = %s and condition_id = %s', [cluster_id, cond_id])
    result = cur.fetchall()
    if len(result) == 0:
        cur.execute('insert into networks_bicluster_conditions (bicluster_id, condition_id) values (%s, %s)',
                    [cluster_id, cond_id])
        return True
    else:
        return False

def add_bicluster_conditions(pgconn, sqliteconn, bicl_map, cond_map, iteration):
    """copy bicluster conditions from sqlite to postgres if not there"""
    print "copying bicluster-conditions memberships..."
    clusters = sorted(bicl_map.keys())
    sqlcur = sqliteconn.cursor()
    for cluster in clusters:
        sqlcur.execute('select name from column_members m join column_names n on m.order_num = n.order_num where iteration = ? and cluster = ?', [iteration, cluster])
        cluster_id = bicl_map[cluster]
        conds = [row[0] for row in sqlcur.fetchall()]
        cond_ids = [cond_map[cond] for cond in conds]
        added = 0
        for cond_id in cond_ids:
            added += 1 if add_bicluster_condition(pgconn, cluster_id, cond_id) else 0
        print "conditions added for cluster ", cluster, ": ", added


def add_motif(pgconn, cluster_id, motifnum, num_sites, evalue, pssm):
    print "add motif for cluster: ", cluster_id, " # sites: ", num_sites
    cur = pgconn.cursor()
    cur.execute('insert into networks_motif (bicluster_id, position, sites, e_value) values (%s, %s, %s, %s) returning id', [cluster_id, motifnum, num_sites, evalue])
    motif_id = cur.fetchone()[0]
    for i, pssm_row in enumerate(pssm):
        # positions are 1-based
        cur.execute('insert into pssms (motif_id, position, a, c, g, t) values (%s, %s, %s, %s, %s, %s)', [motif_id, i + 1, pssm_row[0], pssm_row[1], pssm_row[2], pssm_row[3]])


def add_motifs(pgconn, sqliteconn, bicl_map, iteration):
    """copy bicluster motifs and their pssms"""
    print "copying bicluster motifs..."
    clusters = sorted(bicl_map.keys())
    sqlcur = sqliteconn.cursor()
    cur = pgconn.cursor()

    for cluster in clusters:
        sqlcur.execute('select rowid, motif_num, evalue from motif_infos where iteration = ? and cluster = ?', [iteration, cluster])
        results = sqlcur.fetchall()

        cur.execute('select count(*) from networks_motif where bicluster_id = %s',
                    [bicl_map[cluster]])
        if cur.fetchone()[0] == 0:
            for motif_id, num, evalue in results:
                sqlcur.execute('select count(*) from (select distinct gene_num from motif_annotations where motif_info_id = ?)', [motif_id])
                num_sites = sqlcur.fetchone()[0]
                # extract the pssm
                sqlcur.execute('select a,c,g,t from motif_pssm_rows where motif_info_id = ?',
                               [motif_id])
                pssm = [row for row in sqlcur.fetchall()]
                #print "MOTIF: ", (motif_id, num, evalue)
                #print pssm
                add_motif(pgconn, bicl_map[cluster], num, num_sites, evalue, pssm)
        

if __name__ == '__main__':
    description = 'addnwportal.py - adding a cMonkey/python run to the database'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--rundb', required=True, help='cMonkey output db file')
    parser.add_argument('--ratios', required=True, help='cMonkey ratios file')
    args = parser.parse_args()

    # read the matrix
    matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
    infile = util.read_dfile(args.ratios, has_header=True, quote='\"')
    ratios = matrix_factory.create_from(infile)

    # access the run information
    conn = sqlite3.connect(args.rundb)
    cursor = conn.cursor()
    cursor.execute('select organism, species, num_iterations, num_clusters from run_infos')
    orgcode, species, num_iterations, num_clusters = cursor.fetchone()
    print "organism: %s species: %s iterations: %d clusters: %d" % (orgcode, species,
                                                                    num_iterations,
                                                                    num_clusters)

    # start populating the database
    microbedb, organism = make_microbe(orgcode)
    ncbi_code = microbedb.rsat_info.taxonomy_id
    ucsc_code = UCSC_MAP[orgcode]

    pgconn = psycopg2.connect("dbname='network_portal' user='dj_ango' host='localhost' password='django'")
    pgconn.set_isolation_level(0)  # set auto commit
    species_id = add_species(pgconn, orgcode, species, ncbi_code, ucsc_code)
    network_id = add_network(pgconn, species_id, species)
    biclusters = add_biclusters(pgconn, conn, network_id, num_iterations)
    bicl_map = {bicl[1]: bicl[0] for bicl in biclusters}
    chr_map, genes = add_rsat_genes(pgconn, species_id, microbedb)
    gene_map = {gene[1]: gene[0] for gene in genes}
    thesaurus = organism.thesaurus()
    add_genes = add_ratio_genes(pgconn, species_id, ratios, thesaurus, chr_map, gene_map)
    # rebuild gene map if there were unmapped genes
    if len(add_genes) > 0:
        print "%d genes not in RSAT found in the expression" % len(add_genes)
        genes.extend(add_genes)
        gene_map = {gene[1]: gene[0] for gene in genes}

    conditions = add_conditions(pgconn, network_id, ratios)
    cond_map = {cond[1]: cond[0] for cond in conditions}

    # now we have everything to build bicluster memberships
    add_bicluster_genes(pgconn, conn, thesaurus, bicl_map, gene_map, num_iterations)
    add_bicluster_conditions(pgconn, conn, bicl_map, cond_map, num_iterations)
    add_motifs(pgconn, conn, bicl_map, num_iterations)

    print 'done.'
