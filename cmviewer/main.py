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

def runinfo_factory(cursor, row):
    return RunInfo(*row)

def clusterstat_factory(cursor, row):
    return ClusterStat(*row)

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
        cursor.execute('select median_residual from iteration_stats')
        js_mean_residuals = json.dumps([row[0] for row in cursor.fetchall()])
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

        conn.close()
        cursor= None
        conn = None

        tmpl = env.get_template('index.html')
        progress = "%.2f" % min((runinfo.last_iter / runinfo.num_iters * 100.0), 100.0)
        current_iter = int(iteration)
        return tmpl.render(locals())

def setup_routes():
    d = cherrypy.dispatch.RoutesDispatcher()
    d.connect('main', '/:iteration', controller=ClusterViewerApp(), action="iteration")
    return d

if __name__ == '__main__':
    conf = {'/': {'request.dispatch': setup_routes()},
            '/static': {'tools.staticdir.on': True,
                        'tools.staticdir.dir': os.path.join(current_dir, 'static')}}
    cherrypy.config.update(conf)
    app = cherrypy.tree.mount(None, config=conf)
    cherrypy.quickstart(app)
