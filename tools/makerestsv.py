#!/usr/bin/python
import sqlite3
import argparse
import os.path

"""
This is a script taking a result directory from a run that was run with
--keepmemeout True
and produces tsv files for validation
"""
def meme_to_str(outpath, iteration, cluster):
    lines = ""
    try:
        with open(os.path.join(outpath, 'meme-out-%d-%d' % (iteration, cluster))) as infile:
            lines = [line.replace('\n', '<<<<>>>>') for line in infile.readlines()]
    except:
        pass
    return ''.join(lines)

def write_iteration(conn, outfile, iteration):
    cursor = conn.cursor()
    outfile.write('"cols"\t"dens_string"\t"k"\t"meanp_meme"\t"meme_out"\t"resid"\t"rows"\n')
    for cluster in range(num_clusters):
        cursor.execute('select name from column_members m join column_names c on m.order_num = c.order_num where m.cluster = ? and iteration = ?',
                       [(cluster + 1), iteration])
        colnames = [row[0] for row in cursor.fetchall()]
        cols_out = ",".join(colnames)
        cursor.execute('select score from network_stats where network = \'STRING\' and iteration = ?',
                       [iteration])
        row = cursor.fetchone()
        string_dens = row[0] if row != None else 1.0

        cursor.execute('select pval from motif_stats where seqtype = \'upstream\' and iteration = ?',
                       [iteration])
        row = cursor.fetchone()
        meme_pval = row[0] if row != None else 1.0

        meme_out = meme_to_str(args.cmonkeyout, iteration, cluster + 1)

        cursor.execute('select residual from cluster_residuals where cluster = ? and iteration = ?',
                       [cluster, iteration])
        row = cursor.fetchone()
        resid = row[0] if row != None else 1.0

        cursor.execute('select name from row_members m join row_names r on m.order_num = r.order_num where m.cluster = ? and iteration = ?',
                       [(cluster + 1), iteration])
        rownames = [row[0] for row in cursor.fetchall()]
        rows_out = ",".join(rownames)
        
        outfile.write('"%s"\t%f\t%d\t%f\t"%s"\t%f\t"%s"\n' % (cols_out, string_dens, cluster + 1, meme_pval, meme_out, resid, rows_out))

if __name__ == '__main__':
    description = 'Result maker'
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--cmonkeyout', required=True, help='cmonkey output directory')
    parser.add_argument('--iteration', required=False, default=None,
                        help='iteration number')
    args = parser.parse_args()

    dbfile = '%s/cmonkey_run.db' % args.cmonkeyout

    conn = sqlite3.connect(dbfile)
    cursor = conn.cursor()
    cursor.execute('select num_clusters, last_iteration from run_infos')
    num_clusters, last_iteration = cursor.fetchone()

    if args.iteration != None:
        iteration = int(args.iteration)
        with open("cmresults-%d.tsv" % iteration, 'w') as outfile:
            write_iteration(conn, outfile, iteration)
    else:
        cursor.execute('select iteration from iteration_stats order by iteration')
        iterations = [row[0] for row in cursor.fetchall()]
        for iteration in iterations:
            print "processing iteration ", iteration
            with open("cmresults-%d.tsv" % iteration, 'w') as outfile:
                write_iteration(conn, outfile, iteration)

    conn.close()
