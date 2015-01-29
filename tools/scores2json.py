#!/usr/bin/python
# make sure cmonkey is in PYTHONPATH
import cPickle
import argparse
import json
import sqlite3

def compute_netscores(netscore_file, num_clusters):
    fname = netscore_file + ".json"
    with open(netscore_file) as infile:
        netscores = cPickle.load(infile)

    # compute average
    outscores = {}
    for cluster in range(1, num_clusters + 1):
        score_total = 0.0
        for nwname in netscores.keys():  # network types
            clscore = 0.0
            if cluster in netscores[nwname]:
                for key, score in netscores[nwname][cluster].items():
                    clscore += score
                if clscore != 0.0:
                    clscore /= len(netscores[nwname][cluster].keys())
                score_total += clscore
        outscores[cluster] = score_total / len(netscores.keys())

    with open(fname, "w+") as outfile:
        json.dump(outscores, outfile)


def compute_combscores(combscore_file, num_clusters, row_members):
    fname = combscore_file + ".json"
    with open(combscore_file) as infile:
        combscores = cPickle.load(infile)
    combrows = { row: index for index, row in enumerate(combscores.row_names) }

    # combscores is a DataMatrix object, just extract the relevant scores
    outscores = {}
    for cluster in range(1, num_clusters + 1):
        cluster_scores = { row : combscores.values[combrows[row]][cluster - 1]
                           for row in row_members[cluster] }
        outscores[cluster] = sum(cluster_scores.values()) / len(cluster_scores)

    with open(fname, "w+") as outfile:
        json.dump(outscores, outfile)


if __name__ == '__main__':
    description = """
scores2json - convert cmonkey2 pickle scores to json scores
"""
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--dbfile', required=True,
                        help='database file')
    parser.add_argument('--netscores')
    parser.add_argument('--combscores')
    args = parser.parse_args()

    conn = sqlite3.connect(args.dbfile)
    c = conn.cursor()
    c.execute("""select last_iteration, num_clusters from run_infos""")
    last_iter, num_clusters = c.fetchone()
    print "last iteration: ", last_iter
    print "# clusters: ", num_clusters
    row_members = {}
    for cluster in range(1, num_clusters + 1):
        query = """select rn.name from row_members rm join row_names rn on rm.order_num = rn.order_num where iteration = ? and cluster = ?"""
        c.execute(query, [last_iter, cluster])
        row_members[cluster] = [row[0] for row in c.fetchall()]

    if args.netscores:
        compute_netscores(args.netscores, num_clusters)

    if args.combscores:
        compute_combscores(args.combscores, num_clusters, row_members)
