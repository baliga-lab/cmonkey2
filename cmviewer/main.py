#!/usr/bin/python

import cherrypy
from jinja2 import Environment, FileSystemLoader
import os
import sqlite3
from collections import namedtuple, defaultdict
import json
import gzip
import numpy as np
import glob
import math
import argparse


current_dir = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(current_dir, 'templates')))
outdir = 'out'  # make it flexible

parser = argparse.ArgumentParser()
parser.add_argument('--out', default=outdir, help='output directory')
args = parser.parse_args()
outdir = os.path.join(os.getcwd(), args.out)
outdb = os.path.join(outdir, 'cmonkey_run.db')


RunInfo = namedtuple('RunInfo',
                     ['species', 'orgcode', 'num_iters', 'last_iter',
                      'num_rows', 'num_cols', 'num_clusters', 'start_time', 'finish_time'])

ClusterStat = namedtuple('ClusterStat',
                         ['iter', 'cluster', 'num_rows', 'num_cols', 'residual'])

IterationStat = namedtuple('IterationStat', ['iter', 'label', 'score'])

MotifInfo = namedtuple('MotifInfo', ['id', 'cluster', 'seqtype', 'num', 'evalue'])

MotifPSSMRow = namedtuple('MotifPSSMRow', ['motif_id', 'row', 'a', 'c', 'g', 't'])

MotifAnnotation = namedtuple('MotifAnnotation', ['motif_info_id', 'seqtype', 'motif_num',
                                                 'gene', 'pos', 'reverse', 'pvalue'])

def normalize_js(value):
    if math.isnan(value):
        return 0.0
    else:
        return value
    

def format_float(value):
    if value <= 1e-3:
        return "%.2e" % value
    else:
        return "%.2f" % value

class Ratios:
    """A helper class that provides useful functionality to generate plotting
    related information"""

    def __init__(self, genes, conds, data):
        self.genes = genes
        self.conds = conds
        self.data = data
        self.gene_idx = {gene: i for i, gene in enumerate(genes)}
        self.cond_idx = {cond: i for i, cond in enumerate(conds)}

    def mean(self):
        return np.mean(self.data)


    def subratios_for(self, genes, conds):
        """Arrange cluster expression data for plotting.
        The result is a sub matrix with |genes| rows and the original number of
        columns, but rearranged so that the columns inside the cluster are
        first then the ones that are outside the cluster are last"""
        in_indexes = [self.cond_idx[cond] for cond in conds]
        out_indexes = sorted(set(self.cond_idx.values()) - set(in_indexes))
        col_indexes = in_indexes + out_indexes
        row_indexes = [self.gene_idx[gene] for gene in genes]
        data = self.data[row_indexes][:, col_indexes]
        new_conds = [self.conds[i] for i in col_indexes]
        return Ratios(genes, new_conds, data)

    def hs_subratios_for(self, genes, conds):
        subratios = self.subratios_for(genes, conds)
        return [{'name': gene, 'data': [normalize_js(val) for val in subratios.data[i]]}
                for i, gene in enumerate(genes)]

    def hs_boxplot_data_for(self, genes, conds):
        def make_row(row):
            r = sorted(row)
            minval = r[0]
            maxval = r[-1]
            median = r[len(r) / 2 + len(r) % 2]
            quart = len(r) / 4
            lower_quartile = r[quart]
            upper_quartile = r[-(quart + 1)]
            return [normalize_js(minval), normalize_js(lower_quartile),
                    normalize_js(median), normalize_js(upper_quartile),
                    normalize_js(maxval)]
            
        subratios = self.subratios_for(genes, conds)
        # cut up the data into left and right half
        data_in = subratios.data[:,:len(conds)].T
        data_out = subratios.data[:,len(conds):].T
        inrows = sorted(map(make_row, data_in), key=lambda r: r[2])
        outrows = sorted(map(make_row, data_out), key=lambda r: r[2])
        result = inrows + outrows
        return json.dumps(result)
                

def read_ratios():
    def to_float(s):
        if s == 'NA':
            return float('nan')
        else:
            return float(s)

    ratios_file = os.path.join(outdir, 'ratios.tsv.gz')
    with gzip.open(ratios_file) as infile:
        column_titles = infile.readline().strip().split('\t')
        row_titles = []
        data = []
        for line in infile:
            row = line.strip().split('\t')
            row_titles.append(row[0])
            data.append(map(to_float, row[1:]))
    return Ratios(row_titles, column_titles, np.array(data))

def read_runlogs():
    def read_runlog(fname):
        with open(fname) as infile:
            entries = map(lambda l: 0.0 if l[1] == '1' else float(l[2]),
                          [line.strip().split(':') for line in infile])
        return {'name': os.path.basename(fname).replace('.runlog', ''), 'data': entries}

    return json.dumps([read_runlog(fname)
                       for fname in glob.glob(os.path.join(outdir, '*.runlog'))
                       if os.path.basename(fname) != 'row_scoring.runlog'])

def runinfo_factory(cursor, row):
    return RunInfo(*row)

def motifinfo_factory(cursor, row):
    return MotifInfo(*row)

def motifpssmrow_factory(cursor, row):
    return MotifPSSMRow(*row)

def motifannot_factory(cursor, row):
    return MotifAnnotation(*row)


def clusterstat_factory(cursor, row):
    return ClusterStat(*row)


def iterationstat_factory(cursor, row):
    return IterationStat(*row)


def dbconn():
    return sqlite3.connect(outdb)


def make_int_histogram(counts):
    """input: list of counts
    output: xvalues (count), yvalues (# clusters)"""
    histogram = defaultdict(int)
    for count in counts:
        histogram[count] += 1
    sorted_keys = sorted(histogram.keys())
    return json.dumps(sorted_keys), json.dumps([histogram[key] for key in sorted_keys])
    

def make_float_histogram(values, nbuckets=20):
    minval = min(values)
    maxval = max(values)
    interval = (maxval - minval) / nbuckets
    xvals = ["%.2f" % (minval + interval * i) for i in range(nbuckets)]
    yvals = [0.0] * nbuckets
    for value in values:
        bucketnum = min(nbuckets - 1, int((value - minval) / interval))
        yvals[bucketnum] += 1
    return json.dumps(xvals), json.dumps(yvals)


def make_series(stats):
    groups = defaultdict(list)
    for stat in stats:
        groups[stat.label].append(stat.score)
    return json.dumps([{'name': label, 'data': groups[label]} for label in groups])


class ClusterViewerApp:

    def __init__(self):
        self.__ratios = None

    def ratios(self):
        if self.__ratios is None:
            self.__ratios = read_ratios()
        return self.__ratios

    @cherrypy.expose
    def index(self):
        conn = dbconn()
        cursor = conn.cursor()
        iteration = None
        try:
            cursor.execute('select last_iteration from run_infos')
            iteration = cursor.fetchone()[0]
        except:
            tmpl = env.get_template('not_available.html')
            return tmpl.render(locals())
        finally:
            cursor.close()
            conn.close()
        if iteration is not None:
            raise cherrypy.HTTPRedirect('/%d' % iteration)
        else:
            tmpl = env.get_template('not_available.html')
            return tmpl.render(locals())

    @cherrypy.expose
    def iteration(self, iteration):
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute('select distinct iteration from row_members')
        iterations = [row[0] for row in cursor.fetchall()]
        js_iterations = json.dumps(iterations)
        cursor.close()

        cursor = conn.cursor()
        cursor.execute('select median_residual, fuzzy_coeff from iteration_stats')
        res = [(resid, fuzzy) for resid, fuzzy in cursor.fetchall()]
        js_mean_residuals = json.dumps([row[0] for row in res])
        js_fuzzy_coeff = json.dumps([row[1] for row in res])
        cursor.close()

        conn.row_factory = clusterstat_factory
        cursor = conn.cursor()
        cursor.execute('select iteration, cluster, num_rows, num_cols, residual from cluster_stats')
        row_stats = defaultdict(list)
        col_stats = defaultdict(list)
        resid_stats = defaultdict(list)
        for stat in cursor.fetchall():
            row_stats[stat.iter].append(stat.num_rows)
            col_stats[stat.iter].append(stat.num_cols)
            resid_stats[stat.iter].append(stat.residual)
        js_mean_nrow = json.dumps([float(sum(row_stats[iter])) / len(row_stats[iter])
                                   for iter in sorted(row_stats.keys())])
        js_mean_ncol = json.dumps([float(sum(col_stats[iter])) / len(col_stats[iter])
                                   for iter in sorted(col_stats.keys())])
        js_nrows_x, js_nrows_y = make_int_histogram(row_stats[int(iteration)])
        js_ncols_x, js_ncols_y = make_int_histogram(col_stats[int(iteration)])
        js_resids_x, js_resids_y = make_float_histogram(resid_stats[int(iteration)])
        cursor.close()

        conn.row_factory = runinfo_factory
        cursor = conn.cursor()
        cursor.execute('select species, organism, num_iterations, last_iteration, num_rows, num_columns, num_clusters, start_time, finish_time from run_infos')
        runinfo = cursor.fetchone()
        cursor.close()

        conn.row_factory = iterationstat_factory
        cursor = conn.cursor()
        cursor.execute('select iteration, seqtype, pval from motif_stats')
        js_motif_stats = make_series([row for row in cursor.fetchall()])
        cursor.close()

        cursor = conn.cursor()
        cursor.execute('select iteration, network, score from network_stats')
        js_network_stats = make_series([row for row in cursor.fetchall()])
        cursor.close()

        conn.close()
        cursor= None
        conn = None

        js_runlog_series = read_runlogs()

        tmpl = env.get_template('index.html')
        progress = "%.2f" % min((float(runinfo.last_iter) / float(runinfo.num_iters) * 100.0),
                                100.0)
        current_iter = int(iteration)
        return tmpl.render(locals())

    @cherrypy.expose
    def clusters(self, iteration, *args, **kw):
        print kw
        conn = dbconn()
        conn.row_factory = motifinfo_factory
        cursor = conn.cursor()
        cursor.execute("select rowid, cluster, seqtype, motif_num, evalue from motif_infos where iteration = ?", [iteration])
        # grouped by cluster
        motif_infos = defaultdict(list)
        for row in cursor.fetchall():
            motif_infos[row.cluster].append(row)
        cursor.close()

        conn.row_factory = motifpssmrow_factory
        cursor = conn.cursor()
        cursor.execute("select motif_info_id, row, a, c, g, t from motif_pssm_rows where iteration = ?", [iteration])
        # grouped by motif info id
        motif_pssm_rows = defaultdict(list)
        for row in cursor.fetchall():
            motif_pssm_rows[row.motif_id].append([row.a, row.c, row.g, row.t])
        
        cursor.close()

        conn.row_factory = clusterstat_factory
        cursor = conn.cursor()
        cursor.execute("select iteration, cluster, num_rows, num_cols, residual from cluster_stats where iteration = ? order by residual", [iteration])
        cluster_stats = [row for row in cursor.fetchall()]

        rows = [["%d" % (i + 1),
                 "<a class=\"clusterlink\" id=\"%d\"  href=\"#\">%d</a>" % (stat.cluster, stat.cluster),
                 '%d' % stat.num_rows,
                 '%d' % stat.num_cols,
                 format_float(stat.residual),
                 make_motif_string(motif_infos[stat.cluster],
                                   motif_pssm_rows)]
                for i, stat in enumerate(cluster_stats)]

        cursor.close()
        conn.close()
        return json.dumps({'aaData': rows})

    @cherrypy.expose
    def view_cluster(self, **kw):
        cluster = int(kw['cluster'])
        iteration = int(kw['iteration'])
        conn = dbconn()

        cursor = conn.cursor()
        cursor.execute('select species from run_infos')
        species = cursor.fetchone()[0]
        cursor.close()

        cursor = conn.cursor()
        cursor.execute('select name from row_names rn join row_members rm on rn.order_num = rm.order_num where iteration = ? and cluster = ?', (iteration, cluster))
        rows = [row[0] for row in cursor.fetchall()]
        cursor.close()

        cursor = conn.cursor()
        cursor.execute('select name from column_names cn join column_members cm on cn.order_num = cm.order_num where iteration = ? and cluster = ?', (iteration, cluster))
        columns = [row[0] for row in cursor.fetchall()]
        cursor.close()

        js_ratios = json.dumps(self.ratios().hs_subratios_for(rows, columns))

        # extract motif information
        conn.row_factory = motifinfo_factory
        cursor = conn.cursor()
        cursor.execute("select rowid, cluster, seqtype, motif_num, evalue from motif_infos where iteration = ? and cluster = ?", [iteration, cluster])
        # grouped by seqtype
        motif_infos = defaultdict(list)
        for row in cursor.fetchall():
            motif_infos[row.seqtype].append(row)
        cursor.close()
        seqtypes = motif_infos.keys()

        conn.row_factory = motifpssmrow_factory
        cursor = conn.cursor()
        motif_ids = map(lambda i: str(i.id),
                        [mis for mismis in motif_infos.values() for mis in mismis])
        cursor.execute("select motif_info_id, row, a, c, g, t from motif_pssm_rows where motif_info_id in (%s)" % ','.join(motif_ids))
        # grouped by motif info id
        motif_pssm_rows = defaultdict(list)
        for row in cursor.fetchall():
            motif_pssm_rows[row.motif_id].append([row.a, row.c, row.g, row.t])
        js_motif_pssms = {motif_id: json.dumps({'alphabet':['A','C','G','T'],
                                                'values':motif_pssm_rows[motif_id]})
                          for motif_id in motif_pssm_rows}
        cursor.close()

        # annotations
        conn.row_factory = sqlite3.Row
        cursor = conn.cursor()
        cursor.execute('select motif_info_id, count(row) from motif_pssm_rows where iteration = ? group by motif_info_id', [iteration])
        motif_lengths = {row[0]: row[1] for row in cursor.fetchall()}
        cursor.close()

        conn.row_factory = motifannot_factory
        cursor = conn.cursor()
        cursor.execute('select a.motif_info_id, seqtype, motif_num, g.name, position, reverse, pvalue from motif_annotations a join motif_infos i on a.motif_info_id = i.rowid join row_names g on g.order_num = a.gene_num where i.iteration = ? and i.cluster = ?', [iteration, cluster])
        annotations = [row for row in cursor.fetchall()]
        cursor.close()

        st_annots = defaultdict(list)  # group by seqtype
        for annot in annotations:
            st_annots[annot.seqtype].append(annot)
        st_gene_annots = {}
        for seqtype in st_annots:  # group by gene
            gene_annots = defaultdict(list)
            for annot in st_annots[seqtype]:
                gene_annots[annot.gene].append(annot)
            st_gene_annots[seqtype] = gene_annots

        js_annotation_map = {}
        for seqtype, gene_annots in st_gene_annots.items():
            st_annots = []
            for gene, annots in gene_annots.items():
                matches = [{'motif': a.motif_num - 1, 'start': a.pos,
                            'length': motif_lengths[a.motif_info_id],
                            'reverse': a.reverse, 'score': a.pvalue}
                           for a in annots]
                st_annots.append({'gene': gene, 'condition': '', 'log10': 0.17,
                                  'boxColor': '#08f', 'lineColor': '#000',
                                  'matches': matches})
            js_annotation_map[seqtype] = json.dumps(st_annots)

        conn.close()
        ratios_mean = normalize_js(self.ratios().subratios_for(rows, columns).mean())
        js_boxplot_ratios = self.ratios().hs_boxplot_data_for(rows, columns)
        tmpl = env.get_template('cluster.html')
        return tmpl.render(locals())


def consensus(rows):
    alphabet = ['A', 'C', 'G', 'T']
    result = ''
    for row in rows:
        maxcol = 0
        maxval = 0.0
        for col in range(len(row)):
            if row[col] > maxval:
                maxcol = col
                maxval = row[col]
        result += alphabet[maxcol]
    return result

def make_motif_string(motif_infos, motif_pssm_rows):
    result = "("
    byseqtype = defaultdict(list)
    for info in motif_infos:
        byseqtype[info.seqtype].append(info)
    for seqtype in byseqtype:
        reps = ["%s (%s)" % (consensus(motif_pssm_rows[motif.id]), format_float(motif.evalue))
                for motif in byseqtype[seqtype]]
        result += "%s -> [%s]" % (seqtype, ', '.join(reps))
    result += ")"
    return result

def setup_routes():
    d = cherrypy.dispatch.RoutesDispatcher()
    main = ClusterViewerApp()
    d.connect('main', '/', controller=main, action="index")
    d.connect('iteration', '/:iteration', controller=main, action="iteration")
    d.connect('clusters', '/clusters/:iteration', controller=main, action="clusters")
    d.connect('cluster', '/cluster/:iteration/:cluster', controller=main, action="view_cluster")
    return d

if __name__ == '__main__':
    conf = {'/': {'request.dispatch': setup_routes()},
            '/static': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(current_dir, 'static')}}
    cherrypy.config.update(conf)
    app = cherrypy.tree.mount(None, config=conf)
    cherrypy.quickstart(app)
