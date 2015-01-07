#!/usr/bin/python

import cherrypy
from cherrypy import tools
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
parser.add_argument('--port', type=int, default=8080, help='port to listen to web requests')
args = parser.parse_args()
outdir = os.path.join(os.getcwd(), args.out)
outdb = os.path.join(outdir, 'cmonkey_run.db')


RunInfo = namedtuple('RunInfo',
                     ['species', 'orgcode', 'num_iters', 'last_iter',
                      'num_rows', 'num_cols', 'num_clusters', 'start_time', 'finish_time',
                      'run_secs'])

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
    return sqlite3.connect(outdb, timeout=10, isolation_level=None)


def make_int_histogram(counts):
    """input: list of counts
    output: xvalues (count), yvalues (# clusters)"""
    histogram = defaultdict(int)
    for count in counts:
        histogram[count] += 1
    sorted_keys = sorted(histogram.keys())
    return sorted_keys, [histogram[key] for key in sorted_keys]


def make_float_histogram(values, nbuckets=20):
    if len(values) > 0:
        minval = min(values)
        maxval = max(values)
        interval = (maxval - minval) / nbuckets
        xvals = ["%.2f" % (minval + interval * i) for i in range(nbuckets)]
        yvals = [0.0] * nbuckets
        for value in values:
            bucketnum = min(nbuckets - 1, int((value - minval) / interval))
            yvals[bucketnum] += 1
    else:
        xvals = []
        yvals = []
    return xvals, yvals


def make_series(stats):
    """Creates a data series for Highcharts and returns a triple of the
    series data and minimum and maximum score
    The input is a list of IterationStat objects
    """
    groups = defaultdict(list)
    scores = [stat.score for stat in stats]
    if len(scores) > 0:
        minscore = min(scores)
        maxscore = max(scores)
    else:
        minscore = 0.0
        maxscore = 0.0

    for stat in stats:
        groups[stat.label].append(stat.score)
    minscore = math.floor(minscore)
    maxscore = math.ceil(maxscore)
    return [{'name': label, 'data': groups[label]} for label in groups], minscore, maxscore


class ClusterViewerApp:

    def __init__(self):
        self.__ratios = None

    def ratios(self):
        if self.__ratios is None:
            self.__ratios = read_ratios()
        return self.__ratios

    @cherrypy.expose
    def index(self):
        iteration = None
        conn = None
        cursor = None
        try:
            conn = dbconn()
            cursor = conn.cursor()
            cursor.execute('select max(iteration) from row_members')
            iteration = cursor.fetchone()[0]
        except:
            tmpl = env.get_template('not_available.html')
            return tmpl.render(locals())
        finally:
            if cursor is not None:
                cursor.close()
            if conn is not None:
                conn.close()

        if iteration is not None:
            return self.real_index()
        else:
            tmpl = env.get_template('not_available.html')
            return tmpl.render(locals())

    def filter_clusters(self, cursor, iteration, min_residual, max_residual):
        # used clusters
        query_params = [iteration]
        query = "select distinct cluster from cluster_stats where iteration=?"
        if min_residual is not None:
            query += ' and residual >= ?'
            query_params.append(min_residual)
        if max_residual is not None:
            query += ' and residual <= ?'
            query_params.append(max_residual)
        cursor.execute(query, query_params)
        return [row[0] for row in cursor.fetchall()]

    def filter_motifs(self, cursor, iteration, valid_clusters,
                      min_evalue, max_evalue):
        query_params = [iteration]
        query = "select rowid,cluster,motif_num from motif_infos where iteration=?"
        if min_evalue is not None:
            query += ' and evalue >= ?'
            query_params.append(min_evalue)
        if max_evalue is not None:
            query += ' and evalue <= ?'
            query_params.append(max_evalue)
        cursor.execute(query, query_params)
        return [(mid, cluster, motif_num)
                for mid,cluster,motif_num in cursor.fetchall()
                if cluster in valid_clusters]

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def cytoscape_nodes(self, iteration, min_residual=None, max_residual=None,
                        min_evalue=None, max_evalue=None):
        conn = dbconn()
        cursor = conn.cursor()

        clusters = self.filter_clusters(cursor, iteration, min_residual, max_residual)
        valid_clusters = set(clusters)
        clusters_json = [{'classes': 'clusters',
                          'data': {'id': '%d' % cluster, 'name': '%d' % cluster}}
                         for cluster in clusters]

        # used genes
        cursor.execute("select distinct name,cluster from row_names rn join row_members rm on rn.order_num = rm.order_num where iteration=? order by name", [iteration])
        genes = [row[0] for row in cursor.fetchall() if row[1] in valid_clusters]
        genes_json = [{'classes': 'genes',
                       'data': {'id': '%s' % gene, 'name': '%s' % gene}}
                      for gene in genes]

        # motifs
        motifs = self.filter_motifs(cursor, iteration, valid_clusters,
                                    min_evalue, max_evalue)
        cursor.close()
        conn.close()

        motifs_json = [{'classes': 'motifs',
                        'data': {'id': 'm%d' % m[0],
                                 'name': '%d_%d' % (m[1], m[2])}}
                       for m in motifs]
        return {'nodes': clusters_json + genes_json + motifs_json }


    @cherrypy.expose
    @cherrypy.tools.json_out()
    def cytoscape_edges(self, iteration, min_residual=None, max_residual=None,
                        min_evalue=None, max_evalue=None):
        conn = dbconn()
        cursor = conn.cursor()
        clusters = self.filter_clusters(cursor, iteration, min_residual, max_residual)
        valid_clusters = set(clusters)
        motifs = self.filter_motifs(cursor, iteration, valid_clusters,
                                    min_evalue, max_evalue)
        valid_motifs = {m[0] for m in motifs}

        # edges between clusters and genes
        cursor.execute("select name, cluster from row_members rm join row_names rn on rm.order_num = rn.order_num where iteration=?", [iteration])
        edges = [(row[0], row[1]) for row in cursor.fetchall() if row[1] in valid_clusters]

        # add the edges between motifs and clusters
        for motif in motifs:
            edges.append(("m%d" % motif[0], motif[1]))

        # tomtom (motif-motif) edges
        cursor.execute("select motif_info_id1, motif_info_id2 from tomtom_results ttr join motif_infos mi on ttr.motif_info_id1=mi.rowid where motif_info_id1 <> motif_info_id2 and iteration=?", [iteration])
        for mid1, mid2 in cursor.fetchall():
            if mid1 in valid_motifs and mid2 in valid_motifs:
                edges.append(("m%d" % mid1, "m%d" % mid2))
        cursor.close()
        conn.close()

        return {'edges': [{'data': {'source': '%s' % id1, 'target': '%s' % id2} }
                          for id1, id2 in edges]}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def run_status(self):
        conn = dbconn()
        conn.row_factory = runinfo_factory
        cursor = conn.cursor()
        cursor.execute("select species, organism, num_iterations, last_iteration, num_rows, num_columns, num_clusters, start_time, finish_time, (strftime('%s', finish_time) - strftime('%s', start_time)) from run_infos")
        runinfo = cursor.fetchone()

        cursor.close()
        conn.close()
        progress = "%.2f" % min((float(runinfo.last_iter) / float(runinfo.num_iters) * 100.0),
                                100.0)
        result = {'progress': progress, 'finished': False}
        if runinfo.finish_time:
            elapsed_hours = runinfo.run_secs / 3600
            elapsed_mins = (runinfo.run_secs - (elapsed_hours * 3600)) / 60
            result['elapsed_time'] = "(%d hours %d minutes)" % (elapsed_hours, elapsed_mins)
            result['finish_time'] = runinfo.finish_time
            result['finished'] = True
        return result

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def iterations(self):
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute("select distinct iteration from row_members order by iteration")
        result = [row[0] for row in cursor.fetchall()]
        cursor.close()
        conn.close()
        return result

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def mean_residuals(self):
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute("select score from iteration_stats its join statstypes st on its.statstype = st.rowid where st.name = 'median_residual' order by iteration")
        resids = [row[0] for row in cursor.fetchall()]

        cursor.close()
        conn.close()
        return {'min': min(resids), 'max': max(resids), 'values': resids}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def mean_cluster_members(self):
        row_stats = defaultdict(list)
        col_stats = defaultdict(list)
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute('select iteration, cluster, num_rows, num_cols from cluster_stats')
        for iteration, cluster, nrows, ncols in cursor.fetchall():
            row_stats[iteration].append(nrows)
            col_stats[iteration].append(ncols)
        mean_nrow = [float(sum(row_stats[i])) / len(row_stats[i])
                     for i in sorted(row_stats.keys())]
        mean_ncol = [float(sum(col_stats[i])) / len(col_stats[i])
                     for i in sorted(col_stats.keys())]
        cursor.close()
        conn.close()
        return {'meanNumRows': mean_nrow, 'meanNumCols': mean_ncol}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def runlog(self):
        def read_runlog(fname):
            with open(fname) as infile:
                entries = map(lambda l: 0.0 if l[1] == '1' else float(l[2]),
                              [line.strip().split(':') for line in infile])
            return {'name': os.path.basename(fname).replace('.runlog', ''), 'data': entries}

        return [read_runlog(fname)
                for fname in glob.glob(os.path.join(outdir, '*.runlog'))
                if os.path.basename(fname) not in ['row_scoring.runlog',
                                                   'column_scoring.runlog']]

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def fuzzy_coeffs(self):
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute("select score from iteration_stats its join statstypes st on its.statstype = st.rowid where st.name = 'fuzzy_coeff' order by iteration")
        result = [row[0] for row in cursor.fetchall()]
        cursor.close()
        conn.close()
        return result

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def cluster_row_hist(self):
        """Note: this is actually iteration-specific, currently we lock this to the last
        iteration until it becomes an issue"""
        row_stats = defaultdict(list)
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute('select iteration, cluster, num_rows from cluster_stats')
        for iteration, cluster, nrows in cursor.fetchall():
            row_stats[iteration].append(nrows)
        nrows_x, nrows_y = make_int_histogram(row_stats[max(row_stats.keys())])
        cursor.close()
        conn.close()
        return {'xvalues': nrows_x, 'yvalues': nrows_y}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def cluster_col_hist(self):
        """Note: this is actually iteration-specific, currently we lock this to the last
        iteration until it becomes an issue"""
        col_stats = defaultdict(list)
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute('select iteration, cluster, num_cols from cluster_stats')
        for iteration, cluster, ncols in cursor.fetchall():
            col_stats[iteration].append(ncols)
        ncols_x, ncols_y = make_int_histogram(col_stats[max(col_stats.keys())])
        cursor.close()
        conn.close()
        return {'xvalues': ncols_x, 'yvalues': ncols_y}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def cluster_residuals(self):
        """Note: this is actually iteration-specific, currently we lock this to the last
        iteration until it becomes an issue"""
        resid_stats = defaultdict(list)
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute('select iteration, cluster, residual from cluster_stats')
        for i, cluster, residual in cursor.fetchall():
            resid_stats[i].append(residual)
        resids_x, resids_y = make_float_histogram(resid_stats[max(resid_stats.keys())])
        cursor.close()
        conn.close()
        return {'xvalues': resids_x, 'yvalues': resids_y}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def network_score_means(self):
        conn = dbconn()
        conn.row_factory = iterationstat_factory
        cursor = conn.cursor()
        cursor.execute("select iteration,name,score from iteration_stats its join statstypes st on its.statstype=st.rowid where category='network'")
        series, min_score, max_score = make_series([row for row in cursor.fetchall()])
        cursor.close()
        conn.close()
        return {'min': min_score, 'max': max_score, 'series': series}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def slider_ranges(self, iteration):
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute("select min(residual), max(residual) from cluster_stats where iteration=?", [iteration])
        min_residual, max_residual = cursor.fetchone()
        residual_step = (max_residual - min_residual) / 100.0
        cursor.execute("select min(evalue), max(evalue) from motif_infos where iteration=?",
                       [iteration])
        min_evalue, max_evalue = cursor.fetchone()
        evalue_step = (max_evalue - min_evalue) / 100.0
        cursor.close()
        return {'residual': {'min': min_residual, 'max': max_residual, 'step': residual_step},
                'evalue': {'min': min_evalue, 'max': max_evalue, 'step': evalue_step}}

    @cherrypy.expose
    @cherrypy.tools.json_out()
    def generic_score_means(self):
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute("select rowid,name from statstypes where (category='scoring' or category='seqtype') and name not in ('Rows', 'Columns', 'Networks')")
        types = [row[1] for row in cursor.fetchall()]
        cursor.close()

        conn.row_factory = iterationstat_factory
        cursor = conn.cursor()
        stats_scores = []
        stats = []
        for statstype in types:
            cursor.execute("select iteration, name, score from iteration_stats its join statstypes st on its.statstype = st.rowid where st.name=?", [statstype])
            mean_stats, min_score, max_score = make_series([row for row in cursor.fetchall()])
            stats_scores.append(min_score)
            stats_scores.append(max_score)
            stats.extend(mean_stats)
        min_stats_score = min(stats_scores)
        max_stats_score = max(stats_scores)
        cursor.close()
        conn.close()
        return {'min': min_stats_score, 'max': max_stats_score, 'series': stats}

    def real_index(self):
        conn = dbconn()

        # extract general information about the run
        conn.row_factory = runinfo_factory
        cursor = conn.cursor()
        cursor.execute("select species, organism, num_iterations, last_iteration, num_rows, num_columns, num_clusters, start_time, finish_time, (strftime('%s', finish_time) - strftime('%s', start_time)) from run_infos")
        runinfo = cursor.fetchone()

        cursor.close()
        conn.close()

        tmpl = env.get_template('index.html')
        return tmpl.render(locals())

    @cherrypy.expose
    def clusters(self, iteration, *args, **kw):
        def min_evalue(motif_infos):
            """returns the minimum e-value of the given motif infos"""
            if len(motif_infos) == 0:
                return 0
            else:
                return min([m.evalue for m in motif_infos])

        sort_col = int(kw['iSortCol_0'])
        sort_reverse = kw['sSortDir_0'] == 'desc'
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
        if sort_col == 1:
            cluster_stats = sorted(cluster_stats, key=lambda item: item.cluster,
                                   reverse=sort_reverse)
        elif sort_col == 2:
            cluster_stats = sorted(cluster_stats, key=lambda item: item.num_rows,
                                   reverse=sort_reverse)
        elif sort_col == 3:
            cluster_stats = sorted(cluster_stats, key=lambda item: item.num_cols,
                                   reverse=sort_reverse)
        elif sort_col == 4:
            cluster_stats = sorted(cluster_stats, key=lambda item: item.residual,
                                   reverse=sort_reverse)
        elif sort_col == 5:
            cluster_stats = sorted(cluster_stats, key=lambda item: min_evalue(motif_infos[item.cluster]),
                                   reverse=sort_reverse)
            

        rows = [["%d" % (i + 1),
                 "<a class=\"clusterlink\" id=\"%d\"  href=\"javascript:void(0)\">%d</a>" % (stat.cluster, stat.cluster),
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
        motif_ids = sorted(js_motif_pssms.keys())

        if len(motif_ids) > 0:
            motif1_length = len(motif_pssm_rows[motif_ids[0]])
            motif1_pssm_tsv = "A\tC\tG\tT\n"
            for a, c, g, t in motif_pssm_rows[motif_ids[0]]:
                motif1_pssm_tsv += "%.4f\t%.4f\t%.4f\t%.4f\n" % (a, c, g, t)
        if len(motif_ids) > 1:
            motif2_length = len(motif_pssm_rows[motif_ids[1]])
            motif2_pssm_tsv = "A\tC\tG\tT\n"
            for a, c, g, t in motif_pssm_rows[motif_ids[1]]:
                motif2_pssm_tsv += "%.4f\t%.4f\t%.4f\t%.4f\n" % (a, c, g, t)

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
    if len(motif_infos) == 0:
        return ""
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

    # run status
    d.connect('run_status', '/run_status', controller=main, action="run_status")
    d.connect('iterations', '/iterations', controller=main, action="iterations")

    # highcharts graph value routes
    d.connect('mean_residuals', '/mean_residuals', controller=main, action="mean_residuals")
    d.connect('mean_cluster_members',
              '/mean_cluster_members', controller=main, action="mean_cluster_members")
    d.connect('runlog', '/runlog', controller=main, action="runlog")
    d.connect('fuzzy_coeffs', '/fuzzy_coeffs', controller=main, action="fuzzy_coeffs")
    d.connect('cluster_row_hist', '/cluster_row_hist', controller=main,
              action="cluster_row_hist")
    d.connect('cluster_col_hist', '/cluster_col_hist', controller=main,
              action="cluster_col_hist")
    d.connect('cluster_residuals', '/cluster_residuals', controller=main,
              action="cluster_residuals")
    d.connect('network_score_means', '/network_score_means', controller=main,
              action="network_score_means")
    d.connect('generic_score_means', '/generic_score_means', controller=main,
              action="generic_score_means")
    d.connect('slider_ranges', '/slider_ranges/:iteration', controller=main,
              action="slider_ranges")

    # cytoscape.js routes
    d.connect('cytonodes', '/cytoscape_nodes/:iteration', controller=main,
              action="cytoscape_nodes")
    d.connect('cytoedges', '/cytoscape_edges/:iteration', controller=main,
              action="cytoscape_edges")

    # cluster list and details
    d.connect('clusters', '/clusters/:iteration', controller=main, action="clusters")
    d.connect('cluster', '/cluster/:iteration/:cluster', controller=main, action="view_cluster")
    return d

if __name__ == '__main__':
    conf = {'/': {'request.dispatch': setup_routes()},
            '/static': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(current_dir, 'static')}}
    cherrypy.config.update(conf)
    app = cherrypy.tree.mount(None, config=conf)
    cherrypy.server.socket_host = '0.0.0.0'
    cherrypy.server.socket_port = args.port
    cherrypy.quickstart(app)
