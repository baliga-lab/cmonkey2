#!/usr/bin/python

import cherrypy
from jinja2 import Environment, FileSystemLoader
import os
import sqlite3
from collections import namedtuple, defaultdict
import json

current_dir = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(current_dir, 'templates')))
outdir = os.path.join(os.path.dirname(current_dir), 'out')  # make it flexible
outdb = os.path.join(outdir, 'cmonkey_run.db')

RunInfo = namedtuple('RunInfo',
                     ['species', 'orgcode', 'num_iters', 'last_iter',
                      'num_rows', 'num_cols', 'num_clusters', 'start_time', 'finish_time'])

ClusterStat = namedtuple('ClusterStat',
                         ['iter', 'cluster', 'num_rows', 'num_cols', 'residual'])

IterationStat = namedtuple('IterationStat', ['iter', 'label', 'score'])

MotifInfo = namedtuple('MotifInfo', ['id', 'cluster', 'seqtype', 'num', 'evalue'])

MotifPSSMRow = namedtuple('MotifPSSMRow', ['motif_id', 'row', 'a', 'c', 'g', 't'])


def runinfo_factory(cursor, row):
    return RunInfo(*row)

def motifinfo_factory(cursor, row):
    return MotifInfo(*row)

def motifpssmrow_factory(cursor, row):
    return MotifPSSMRow(*row)


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

        tmpl = env.get_template('index.html')
        progress = "%.2f" % min((runinfo.last_iter / runinfo.num_iters * 100.0), 100.0)
        current_iter = int(iteration)
        return tmpl.render(locals())

    @cherrypy.expose
    def clusters(self, iteration, *args, **kw):
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
                 '%.2e' % stat.residual,            
                 make_motif_string(motif_infos[stat.cluster],
                                   motif_pssm_rows)]
                for i, stat in enumerate(cluster_stats)]

        cursor.close()
        conn.close()
        return json.dumps({'aaData': rows})


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
        reps = ["%s (%.2e)" % (consensus(motif_pssm_rows[motif.id]), motif.evalue)
                for motif in byseqtype[seqtype]]
        result += "%s -> [%s]" % (seqtype, ', '.join(reps))
    result += ")"
    return result

def setup_routes():
    d = cherrypy.dispatch.RoutesDispatcher()
    main = ClusterViewerApp()
    d.connect('main', '/:iteration', controller=main, action="iteration")
    d.connect('main', '/clusters/:iteration', controller=main, action="clusters")
    return d

if __name__ == '__main__':
    conf = {'/': {'request.dispatch': setup_routes()},
            '/static': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(current_dir, 'static')}}
    cherrypy.config.update(conf)
    app = cherrypy.tree.mount(None, config=conf)
    cherrypy.quickstart(app)

"""
select cluster, residual, num_rows, num_cols from cluster_stats where iteration = 2000;
"""
