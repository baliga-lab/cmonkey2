#!/usr/bin/python
import sys
import csv
import argparse
import sqlite3


def read_regulondb(regulondb_filename, tfs, strong_evidence):
    # TF name, but use common name for regulated gene, because of ratios matrix
    with open(regulondb_filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
        regulondb_rows = [(row[0], row[4], row[7]) for row in reader]

    regulondb_rows = regulondb_rows[1:]
    num_before = len(regulondb_rows)
    regulondb_rows = [row for row in regulondb_rows if row[0] in tfs]
    if strong_evidence:
        regulondb_rows = [row for row in regulondb_rows if row[2] == 'TRUE']

    num_after = len(regulondb_rows)
    #print "# RegulonDB rows before: %d # after: %d" % (num_before, num_after)
    return regulondb_rows


def read_tomtom(tomtom_filename, pval_threshold):
    with open(tomtom_filename) as csvfile:
        reader = csv.reader(csvfile, delimiter='\t', quotechar='\"')
        tomtom_rows = [(row[0], row[1], row[3], row[4]) for row in reader]
    num_before = len(tomtom_rows) - 1    
    tomtom_rows = [(row[0], row[1], float(row[2]), float(row[3]))
                   for row in tomtom_rows[1:] if float(row[2]) < pval_threshold]
    num_after = len(tomtom_rows)
    #print "# TOMTOM results before: %d # after: %d" % (num_before, num_after)
    return tomtom_rows


def make_tfs2motifs(tomtom_rows):
    gene2motifs = {}
    for row in tomtom_rows:
        motif = row[0]
        gene = row[1]
        if not gene in gene2motifs:
            gene2motifs[gene] = set()
        gene2motifs[gene].add(motif)
    return gene2motifs

def genes_for_motifs(conn, motifs, gene_map):
    cursor = conn.cursor()
    result = []
    for motif in motifs:
        comps = motif.split('_')
        cluster = comps[1]
        motifnum = comps[2]
        cursor.execute('select rowid from motif_infos where iteration = 2001 and cluster = ? and motif_num = ?', [cluster, motifnum])
        motif_id = cursor.fetchone()[0]
        cursor.execute('select distinct gene_num from motif_annotations where motif_info_id = ?',
                       [motif_id])
        genes = [gene_map[int(row[0])] for row in cursor.fetchall()]
        result.extend(genes)
    cursor.close()
    return set(result)

def make_tfs2genes(tfs2motifs, cmonkeyout):
    result = {}
    conn = sqlite3.connect(cmonkeyout)
    cursor = conn.cursor()
    cursor.execute('select order_num, name from row_names')
    gene_map = {int(row[0]): row[1] for row in cursor.fetchall()}
    cursor.close()
    for tfs, motifs in tfs2motifs.items():
        result[tfs] = genes_for_motifs(conn, motifs, gene_map)
    conn.close()
    return result

def compute_precision_recall(regulondb_filename, tomtom_filename, pval_threshold):
    # real work starts here
    tomtom_rows = read_tomtom(tomtom_filename, pval_threshold)
    tfs2motifs = make_tfs2motifs(tomtom_rows)
    tfs = list(tfs2motifs.keys())

    regulondb_rows = read_regulondb(regulondb_filename, tfs, args.strongevidence)
    regulondb_tfs2genes = {}
    for row in regulondb_rows:
        tf = row[0]
        gene = row[1]
        if tf not in regulondb_tfs2genes:
            regulondb_tfs2genes[tf] = set()
        regulondb_tfs2genes[tf].add(gene)

    #print regulondb_tfs2genes
    tfs2genes = make_tfs2genes(tfs2motifs, args.cmonkeyout)
    size_before = sum([len(tfs2genes[tf]) for tf in tfs2genes])
    tfs2genes = {tf:tfs2genes[tf] for tf in tfs2genes if tf in regulondb_tfs2genes}
    size_after = sum([len(tfs2genes[tf]) for tf in tfs2genes])
    #print "# tfs before: %d after: %d" % (size_before, size_after)

    total_regulondb = sum([len(regulondb_tfs2genes[tf])
                           for tf in regulondb_tfs2genes])
    true_positives = 0
    false_positives = 0

    for tf, check_genes in tfs2genes.items():
        rdb_genes = regulondb_tfs2genes[tf]
        for gene in check_genes:
            if gene in rdb_genes:
                true_positives += 1
            else:
                false_positives += 1

    total_predictions = true_positives + false_positives
    #print "TP = %d FP = %d, total predictions = %d" % (true_positives, false_positives,
    #                                                   total_predictions)
    precision = float(true_positives) / float(total_predictions)
    recall = float(true_positives) / float(total_regulondb)
    return (precision, recall)

if __name__ == '__main__':
    description = 'Precision/Recall cMonkey/Python'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--regulondb', required=True, help='RegulonDB CSV')
    parser.add_argument('--tomtom', required=True, help='TomTom results')
    parser.add_argument('--cmonkeyout', required=True, help='cMonkey sqlite output')
    parser.add_argument('--threshold', help='p-value threshold', default='0.1')
    parser.add_argument('--strongevidence', action='store_true',
                        help='use only strong evidence',
                        default=True)
    args = parser.parse_args()

    thresholds = [0.1 / (step * 10) for step in range(1, 20)]
    #threshold = float(args.threshold)
    print "threshold\tprecision\trecall"
    for threshold in thresholds:
        precision, recall = compute_precision_recall(args.regulondb, args.tomtom,
                                                     threshold)
        #print "Precision: %f, Recall: %f" % (precision, recall)
        print "%f\t%f\t%f" % (threshold, precision, recall)

