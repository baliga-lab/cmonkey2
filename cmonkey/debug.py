import re
import os


############################################################
#### DEBUGGING
############################################################


def get_last_meme_iteration(outdir):
    """Returns the last iteration where MEME was run"""
    pat = re.compile('meme-out-(\\d+)-(\\d+)$')
    filenames = os.listdir(outdir)
    iterations = []
    for name in filenames:
        m = pat.match(name)
        if m:
            iterations.append(int(m.group(1)))
    if len(iterations) > 0:
        return max(iterations)
    return None


def meme_to_str(outdir, iteration, cluster):
    """get the meme out file as a string"""
    lines = ""
    try:
        with open(os.path.join(outdir, 'meme-out-%04d-%04d' % (iteration, cluster))) as infile:
            lines = [line.replace('\n', '<<<<>>>>') for line in infile.readlines()]
    except:
        pass
    return ''.join(lines)


def write_iteration(conn, outfile, iteration, num_clusters, outdir):
    """writes the iteration into a debug file"""
    outfile.write('"cols"\t"dens_string"\t"k"\t"meanp_meme"\t"meme_out"\t"resid"\t"rows"\n')
    for cluster in range(1, num_clusters + 1):
        cursor = conn.cursor()
        cursor.execute('select name from column_members m join column_names c on m.order_num = c.order_num where m.cluster = ? and iteration = ?', [cluster, iteration])
        colnames = [row[0] for row in cursor.fetchall()]
        cols_out = ",".join(colnames)
        cursor.close()

        cursor = conn.cursor()
        cursor.execute("select score from iteration_stats its join statstypes st on st.rowid=its.statstype where category='network' and name='STRING' and iteration=?",
                       [iteration])
        row = cursor.fetchone()
        string_dens = row[0] if row != None else 1.0
        cursor.close()

        cursor = conn.cursor()
        cursor.execute("select score from iteration_stats its join statstypes st on st.rowid=its.statstype where category='seqtype' and name='upstream' and iteration=?",
                       [iteration])
        row = cursor.fetchone()
        meme_pval = row[0] if row != None else 1.0

        last_meme_iteration = get_last_meme_iteration(outdir)
        meme_out = meme_to_str(outdir, last_meme_iteration, cluster)
        cursor.close()

        cursor = conn.cursor()
        cursor.execute('select residual from cluster_stats where cluster = ? and iteration = ?',
                       [cluster, iteration])
        row = cursor.fetchone()
        resid = row[0] if row != None else 1.0
        cursor.close()

        cursor = conn.cursor()
        cursor.execute('select name from row_members m join row_names r on m.order_num = r.order_num where m.cluster = ? and iteration = ?', [cluster, iteration])
        rownames = [row[0] for row in cursor.fetchall()]
        rows_out = ",".join(rownames)
        cursor.close()
        
        outfile.write('"%s"\t%f\t%d\t%f\t"%s"\t%f\t"%s"\n' % (cols_out, string_dens, cluster, meme_pval, meme_out, resid, rows_out))
