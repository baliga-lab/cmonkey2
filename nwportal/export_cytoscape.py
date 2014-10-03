#!/usr/bin/python

import networkx as nx
import sqlite3
import argparse
import os.path
import psycopg2


###############################################################################
##### Helpers
###############################################################################

def write_graph(org, graph):
    """dumps a NetworkX graph into a GraphML file"""
    writer = nx.readwrite.graphml.GraphMLWriter(encoding='utf-8',prettyprint=True)
    writer.add_graph_element(graph)
    with open('%s-network.xml' % org, 'w') as outfile:
        writer.dump(outfile)

###############################################################################
##### cMonkey/Python export
###############################################################################

def export_cmonkey_results(resultdir, org):
    """exports a network from a cMonkey/Python result file"""
    resultdb = os.path.join(resultdir, 'cmonkey_run.db')
    conn = sqlite3.connect(resultdb)
    cursor = conn.cursor()
    cursor.execute('select last_iteration from run_infos')
    iteration = cursor.fetchone()[0]
    print "fetching iteration ", iteration
    cursor.execute('select distinct cluster from row_members where iteration = ?',
                   [iteration])
    clusters = sorted([row[0] for row in cursor.fetchall()])

    cursor.execute('select distinct cluster, residual from cluster_residuals where iteration = ?',
                   [iteration])
    cluster_residuals = { cluster: residual for cluster, residual in cursor.fetchall() }

    cursor.execute('select name from row_names order by order_num')
    row_names = [row[0] for row in cursor.fetchall()]

    #print clusters

    graph = nx.Graph()
    used_genes = set()
    cluster_genes = {}
    for cluster in clusters:
        graph.add_node("bicluster:%d" % cluster,
                       {'type': 'bicluster',
                        'name': "Bicluster %d" % cluster,
                        'residual': str(cluster_residuals[cluster])
                        })

        cursor.execute("select order_num from row_members where iteration = ? and cluster = ?",
                       [iteration, cluster])
        genes = sorted([row_names[row[0]] for row in cursor.fetchall()])
        cluster_genes[cluster] = genes
        used_genes.update(genes)
        
        #print "cluster ", cluster, " -> ", len(genes), " genes"
    #print used_genes
    for gene in used_genes:
        graph.add_node(gene, {'type':'gene', 'name': gene })


    cursor.execute('select distinct cluster, motif_num, evalue from motif_infos where iteration = ?',
                   [iteration])
    for cluster, motif_num, evalue in cursor.fetchall():
        motif_id = "motif_%d-%d" % (cluster, motif_num)
        graph.add_node(motif_id, {'type':'motif',
                                  'name': motif_id,
                                  'evalue': str(evalue)
                                  })
        graph.add_edge("bicluster:%d" % cluster, motif_id)


    # edges
    for cluster in clusters:
        for gene in cluster_genes[cluster]:
            graph.add_edge("bicluster:%d" % cluster, gene)

    write_graph(org, graph)

###############################################################################
##### Network portal export
###############################################################################

def nwportal_networkid(conn, org):
    cursor = conn.cursor()
    cursor.execute('select nw.id from networks_species sp join networks_network nw on sp.id = nw.species_id where sp.short_name = %s', [org])
    return cursor.fetchone()[0]

def nwportal_clusters(conn, nwid):
    cursor = conn.cursor()
    cursor.execute('select k, residual from networks_bicluster where network_id = %s', [nwid])
    return [(k, residual) for k, residual in cursor.fetchall()]
    

def export_nwportal_network(conn, org):
    """exports a network from the Network Portal database"""
    cursor = conn.cursor()

    graph = nx.Graph()
    used_genes = set()
    cluster_genes = {}

    nwid = nwportal_networkid(conn, org)
    clusters = nwportal_clusters(conn, nwid)

    for cluster, residual in clusters:
        graph.add_node("bicluster:%d" % cluster,
                       {'type': 'bicluster',
                        'name': "Bicluster %d" % cluster,
                        'residual': str(residual)
                        })

    cursor.execute('select k, name  from networks_bicluster bc join networks_bicluster_genes bg on bc.id = bg.bicluster_id join networks_gene g on bg.gene_id = g.id where network_id = %s order by k', [nwid])

    for cluster, gene in cursor.fetchall():
        used_genes.update(gene)
        if cluster not in cluster_genes:
            cluster_genes[cluster] = []
        cluster_genes[cluster].append(gene)

    for gene in used_genes:
        graph.add_node(gene, {'type':'gene', 'name': gene })

    for cluster, residual in clusters:
        for gene in cluster_genes[cluster]:
            graph.add_edge("bicluster:%d" % cluster, gene)

    cursor.execute('select k, position, e_value from networks_bicluster bc join networks_motif m on bc.id = m.bicluster_id where network_id = %s', [nwid])
    cluster_motifs = {}

    for cluster, motif_num, evalue in cursor.fetchall():
        motif_id = "motif_%d-%d" % (cluster, motif_num)
        graph.add_node(motif_id, {'type':'motif',
                                  'name': motif_id,
                                  'evalue': str(evalue)
                                  })
        graph.add_edge("bicluster:%d" % cluster, motif_id)

    write_graph(org, graph)


def export_nwportal_expression(conn, org):
    nwid = nwportal_networkid(conn, org)
    cursor = conn.cursor()
    gene_query = 'select distinct bg.gene_id from networks_bicluster bc join networks_bicluster_genes bg on bc.id = bg.bicluster_id where network_id = %d' % nwid
    cond_query = 'select distinct bcc.condition_id from networks_bicluster bc join networks_bicluster_conditions bcc on bc.id = bcc.bicluster_id where network_id = %d' % nwid
    expr_query = ('select gene_id, condition_id, value from expression where gene_id in (%s) or condition_id in (%s)' % (gene_query, cond_query))

    gene_name_query = 'select distinct g.id, g.name from networks_bicluster bc join networks_bicluster_genes bg on bc.id = bg.bicluster_id join networks_gene g on g.id = bg.gene_id where network_id = %d' % nwid
    cond_name_query = 'select distinct c.id, c.name from networks_bicluster bc join networks_bicluster_conditions bcc on bc.id = bcc.bicluster_id join networks_condition c on c.id = bcc.condition_id where bc.network_id = %d' % nwid

    cursor.execute(gene_query)
    gene_ids = sorted([row[0] for row in cursor.fetchall()])
    cursor.execute(gene_name_query)
    gene_names = {gene_id: name for gene_id, name in cursor.fetchall()}

    cursor.execute(cond_query)
    cond_ids = sorted([row[0] for row in cursor.fetchall()])
    cursor.execute(cond_name_query)
    cond_names = {cond_id: name for cond_id, name in cursor.fetchall()}
    condition_names = [cond_names[cond_id] for cond_id in cond_ids]

    print "# genes: ", len(gene_ids)
    print "# conds: ", len(cond_ids)

    cursor.execute(expr_query)
    expr_values = {(gene_id, cond_id): value for gene_id, cond_id, value in cursor.fetchall()}
    print "# values: ", len(expr_values), " - should ideally be: ", (len(gene_ids) * len(cond_ids))
    with open('%s-ratios.tsv' % org, 'w') as outfile:
        outfile.write("Locus\t")
        outfile.write("\t".join(condition_names))
        outfile.write("\n")
        for gene_id in gene_ids:
            outfile.write(gene_names[gene_id] + "\t")
            values = []
            for cond_id in cond_ids:
                key = (gene_id, cond_id)
                values.append(str(expr_values.get(key, 'NA')))
            outfile.write('\t'.join(values))
            outfile.write('\n')


if __name__ == '__main__':
    description = 'addnwportal.py - adding a cMonkey/python run to the database'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--basedir', required=False, help='cMonkey result directory')
    args = parser.parse_args()
    if args.basedir != None:
        for org in ['bce', 'bsu', 'bth', 'cac', 'cje', 'gsu', 'pae', 'rsp']:
            resultdir = os.path.join(args.basedir, org)
            print "processing directory: ", resultdir
            export_cmonkey_results(resultdir, org)
    else:
        print "export from network portal"
        conn = psycopg2.connect("dbname=network_portal user=dj_ango password=django")
        for org in ['dvu', 'hal', 'mmp', 'syf']:
            export_nwportal_network(conn, org)
            export_nwportal_expression(conn, org)
