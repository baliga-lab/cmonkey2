#!/usr/bin/python

import cherrypy
from jinja2 import Environment, FileSystemLoader
import os
import sqlite3
from collections import namedtuple
import json

current_dir = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(current_dir, 'templates')))
outdir = os.path.join(os.path.dirname(current_dir), 'out')  # make it flexible
outdb = os.path.join(outdir, 'cmonkey_run.db')

RunInfo = namedtuple('RunInfo',
                     ['species', 'orgcode', 'num_iters', 'last_iter',
                      'num_rows', 'num_cols', 'num_clusters', 'start_time', 'finish_time'])

def runinfo_factory(cursor, row):
    return RunInfo(*row)

def dbconn():
    return sqlite3.connect(outdb)


class ClusterViewerApp:

    @cherrypy.expose
    def iteration(self, iteration):
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute('select distinct iteration from row_members')
        iterations = [row[0] for row in cursor.fetchall()]
        cursor.close()

        cursor = conn.cursor()
        cursor.execute('select median_residual from iteration_stats')
        mean_residuals = [row[0] for row in cursor.fetchall()]
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
        js_iterations = json.dumps(iterations)
        js_mean_residuals = json.dumps(mean_residuals)
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
