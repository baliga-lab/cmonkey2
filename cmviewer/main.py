#!/usr/bin/python

import cherrypy
from jinja2 import Environment, FileSystemLoader
import os
import sqlite3
from collections import namedtuple

current_dir = os.path.dirname(os.path.abspath(__file__))
env = Environment(loader=FileSystemLoader(os.path.join(current_dir, 'templates')))
outdir = os.path.join(os.path.dirname(current_dir), 'out')  # make it flexible
outdb = os.path.join(outdir, 'cmonkey_run.db')

RunInfo = namedtuple('RunInfo',
                     ['species', 'orgcode', 'num_iters', 'last_iter',
                      'num_rows', 'num_cols', 'num_clusters'])
def runinfo_factory(cursor, row):
    return RunInfo(*row)

def dbconn():
    return sqlite3.connect(outdb)


class ClusterViewerApp:

    @cherrypy.expose
    def iteration(self, iteration):
        tmpl = env.get_template('index.html')
        conn = dbconn()
        cursor = conn.cursor()
        cursor.execute('select distinct iteration from row_members')
        iterations = [row[0] for row in cursor.fetchall()]
        cursor.close()
        conn.row_factory = runinfo_factory
        cursor = conn.cursor()
        cursor.execute('select species, organism, num_iterations, last_iteration, num_rows, num_columns, num_clusters from run_infos')
        runinfo = cursor.fetchone()
        cursor.close()
        conn.close()
        return tmpl.render(salutation='Hello', target='World',
                           runinfo=runinfo, iterations=iterations)

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
