#!/usr/bin/python
import argparse
import sqlite3
import cmonkey.organism as org
import cmonkey.util
import cmonkey.rsat
import cmonkey.microbes_online
import psycopg2
import collections
import cmonkey.datamatrix as dm
import os.path
import re

# Script to add a cMonkey python run to the network portal database
KEGG_FILE_PATH = 'testdata/KEGG_taxonomy'
GO_FILE_PATH = 'testdata/proteome2taxid'
RSAT_BASE_URL = 'http://rsat.ccb.sickkids.ca'
COG_WHOG_URL = 'ftp://ftp.ncbi.nih.gov/pub/COG/COG/whog'
STRING_URL_PATTERN = "http://networks.systemsbiology.net/string9/%s.gz"
CACHE_DIR = 'cache'
# Map for UCSC genome browser
UCSC_MAP = {'gsu': 'geobSulf', 'cac': '', 'bce': '', 'bsu': 'baicSubt2',
            'rsp': 'rhodSpha', 'cje': 'campJeju', 'ype': 'yersPest_CO92',
            'bth': 'bactThet_VPI_5482', 'ttj': 'therTher_HB8',
            'pae': 'pseuAeru', 'eco': 'eschColi_K12'
             }
PSQL = "dbname='network_portal' user='dj_ango' host='localhost' password='django'"

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


def replace_underscores(input):
    """a way to get those underscores out of the RSAT name and make them
    compatible to the scheme in NCBI: attempt to match a version numbered
    name, and replace the underscores in the version number with dots,
    otherwise, simply, replace all underscores with spaces"""
    regex = re.compile('([^0-9]+)\d(_\d)*')
    if regex.match(input):
        prefix = regex.group(1)
        return prefix.replace('_', ' ') + input.replace(prefix, '').replace('_', '.')
    else:
        return input.replace('_', ' ')
    

def add_species(pgconn, orgcode, species, ncbi_code, ucsc_code):
    """add a species entry for the specified organism"""
    cur = pgconn.cursor()
    cur.execute('select id from networks_species where short_name = %s', [orgcode])
    result = cur.fetchall()
    if len(result) == 0:
        # TODO: the name should not have underscores, currently, what comes
        # out of the organism mapper has underscores
        print "adding species: ", orgcode
        cur.execute('insert into networks_species (name, short_name, ncbi_taxonomy_id, ucsc_id, created_at) values (%s, %s, %s, %s, now()) returning id',
                    [replace_underscores(species), orgcode, ncbi_code, ucsc_code])
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

CHROMOSOME_LENGTHS = {
    'NC_002939.4': 3814139,
    'NC_002516.2': 6264404,
    'NC_004722.1': 5411809,
    'NC_004721.2': 15274,
    'NC_000964.2': 4214630,
    'NC_004663.1': 6260361,
    'NC_004703.1': 33038,
    'NC_003030.1': 3940880,
    'NC_001988.2': 192000,
    'NC_002163.1': 1641481,
    'NC_007493.1': 3188609,
    'NC_009007.1': 114045,
    'NC_007488.1': 114178,
    'NC_007490.1': 100828,
    'NC_009008.1': 37100,
    'NC_007489.1': 105284,
    'NC_007494.1': 943016,
    'NC_000913.2': 4639675
}

CHROMOSOME_NAMES = {
    'NC_002939.4': 'chromosome',
    'NC_002516.2': 'chromosome',
    'NC_004722.1': 'chromosome',
    'NC_004721.2': 'plasmid',
    'NC_000964.2': 'chromosome',
    'NC_004663.1': 'chromosome',
    'NC_004703.1': 'plasmid',
    'NC_003030.1': 'chromosome',
    'NC_001988.2': 'plasmid',
    'NC_002163.1': 'chromosome',
    'NC_007493.1': 'chromosome 1',
    'NC_009007.1': 'plasmid A',
    'NC_007488.1': 'plasmid B',
    'NC_009008.1': 'plasmid E',
    'NC_007490.1': 'plasmid D',
    'NC_007489.1': 'plasmid C',
    'NC_007494.1': 'chromosome 2',
    'NC_000913.2': 'chromosome'
}

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
    """add the specified chromosome to the database or return its entry
    Note: this does not try to retrieve the chromosome length at the moment,
    we have to do it by hand"""
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

def add_bicluster_genes(pgconn, sqliteconn, species_id, chr_map, thesaurus, bicl_map, gene_map, iteration):
    """copy biclusters from sqlite to postgres if not there"""
    print "copying bicluster-gene memberships..."
    chr_id = chr_map.items()[0][1]
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
                primary_name = thesaurus[gene]
                #print "%s in thesaurus -> %s" % (gene, primary_name)
                mgenes.append(primary_name)
                # handle the weird case that a gene is in the thesaurus
                # but does not have a feature entry: in this case,
                # we need to add a dummy gene and add the entry to the gene map
                if primary_name not in gene_map:
                    print "WARNING: weird RSAT error, need to add a dummy gene %s" % primary_name
                    gene_pair = add_gene(pgconn, species_id, chr_id, primary_name, primary_name,
                                       'DUMMY', 0, 0, '+')
                    gene_map[primary_name] = gene_pair[0]
            else:
                if gene not in gene_map:
                    print "WARNING: weird gene-RSAT inconsistency, need to add a dummy gene %s" % gene
                    gene_pair = add_gene(pgconn, species_id, chr_id, gene, gene,
                                         'DUMMY', 0, 0, '+')
                    gene_map[gene] = gene_pair[0]
                mgenes.append(gene)
        gene_ids = [gene_map[gene] for gene in mgenes]
        added = 0
        for gene_id in gene_ids:
            added += 1 if add_bicluster_gene(pgconn, cluster_id, gene_id) else 0
        if added > 0:
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
        if added > 0:
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


def add_motif_annotations(pgconn, sqliteconn, bicl_map, gene_map, thesaurus, iteration):
    """copy motif annotations"""
    print "copying motif annotations"
    sqlcur = sqliteconn.cursor()
    sqlcur.execute('select name from row_names order by order_num')
    src_genes = [row[0] for row in sqlcur.fetchall()]
    clusters = sorted(bicl_map.keys())
    cur = pgconn.cursor()
    for cluster in clusters:
        sqlcur.execute('select rowid, motif_num from motif_infos where cluster = ? and iteration = ?',
                       [cluster, iteration])
        motif_infos = [row for row in sqlcur.fetchall()]
        for src_motif_id, motif_num in motif_infos:
            cur.execute('select id from networks_motif where bicluster_id = %s and position = %s', [bicl_map[cluster], motif_num])
            dest_motif_id = cur.fetchone()[0]

            cur.execute('select count(*) from networks_motifannotation where motif_id = %s',
                        [dest_motif_id])
            num_dest_annots = cur.fetchone()[0]

            sqlcur.execute('select gene_num, position, reverse, pvalue from motif_annotations where motif_info_id = ?', [src_motif_id])
            annots = [row for row in sqlcur.fetchall()]
            if num_dest_annots == 0 and len(annots) > 0:
                print "add annotations, cluster ", cluster, " motif num: ", motif_num, " src id = ", src_motif_id, " dest id = ", dest_motif_id, " # annotations: ", len(annots)
                for gene_num, position, reverse, pvalue in annots:
                    gene = src_genes[gene_num]
                    if gene in thesaurus:
                        gene = thesaurus[gene]
                    gene_id = gene_map[gene]
                    reverse = True if reverse == 1 else False
                    cur.execute('insert into networks_motifannotation (motif_id, gene_id, position, reverse, pvalue) values (%s,%s,%s,%s,%s)', [dest_motif_id, gene_id, position, reverse, pvalue])
                    pass
            

def add_expressions(pgconn, ratios, thesaurus, gene_map, cond_map, exptable):
    """generate the gene expression table. Gene expressions are huge and
    inserting them separately is very slow. Since we use Postgres, the
    fastest way is to write out a tab-separated file and use the copy
    command to add everything
    """
    cur = pgconn.cursor()
    # select some items to determine whether we have any gene expressions
    genes = ratios.row_names[:5]
    gene_ids = []
    for gene in genes:
        if gene in thesaurus:
            gene_ids.append(gene_map[thesaurus[gene]])
        else:
            gene_ids.append(gene_map[gene])
    
    with open(exptable, 'w') as outfile:        
        for row in range(ratios.num_rows):
            gene = ratios.row_names[row]
            if gene in thesaurus:
                gene_id = gene_map[thesaurus[gene]]
            else:
                gene_id = gene_map[gene]

            for col in range(ratios.num_columns):
                cond_id = cond_map[ratios.column_names[col]]
                value = ratios.values[row][col]
                outfile.write("%d\t%d\t%f\n" % (gene_id, cond_id, value))

if __name__ == '__main__':
    description = 'addnwportal.py - adding a cMonkey/python run to the database'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--resultdir', required=True, help='cMonkey result directory')
    parser.add_argument('--exptable', help='filename of expression table to generate',
                        default=None)
    args = parser.parse_args()
    resultdb = os.path.join(args.resultdir, 'cmonkey_run.db')
    ratiofile = os.path.join(args.resultdir, 'ratios.tsv.gz')

    # read the matrix
    matrix_factory = dm.DataMatrixFactory([dm.nochange_filter, dm.center_scale_filter])
    infile = util.read_dfile(ratiofile, has_header=True, quote='\"')
    ratios = matrix_factory.create_from(infile)

    # access the run information
    conn = sqlite3.connect(resultdb)
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

    pgconn = psycopg2.connect(PSQL)
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
    add_bicluster_genes(pgconn, conn, species_id, chr_map, thesaurus, bicl_map, gene_map, num_iterations)
    add_bicluster_conditions(pgconn, conn, bicl_map, cond_map, num_iterations)
    add_motifs(pgconn, conn, bicl_map, num_iterations)

    if args.exptable != None:
        print "writing expression table..."
        add_expressions(pgconn, ratios, thesaurus, gene_map, cond_map, args.exptable)

    add_motif_annotations(pgconn, conn, bicl_map, gene_map,
                          thesaurus, num_iterations)

    print 'done.'
