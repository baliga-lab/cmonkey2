import re
import os
import cmonkey.database as cm2db
from sqlalchemy import and_


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


def write_iteration(session, outfile, iteration, num_clusters, outdir, as_binary=True):
    """writes the iteration into a debug file"""
    HEADER = '"cols"\t"dens_string"\t"k"\t"meanp_meme"\t"meme_out"\t"resid"\t"rows"\n'
    if as_binary:
        outfile.write(HEADER.encode('utf-8'))
    else:
        outfile.write(HEADER)

    for cluster in range(1, num_clusters + 1):
        colnames = [colmemb.column_name.name
                    for colmemb in session.query(cm2db.ColumnMember).filter(and_(cm2db.ColumnMember.cluster == cluster,
                                                                                 cm2db.ColumnMember.iteration == iteration))]
        cols_out = ','.join(colnames)

        rownames = [rowmemb.row_name.name
                    for rowmemb in session.query(cm2db.RowMember).filter(and_(cm2db.RowMember.cluster == cluster,
                                                                              cm2db.RowMember.iteration == iteration))]
        rows_out = ','.join(rownames)

        try:
            string_dens = session.query(cm2db.IterationStat).join(
                cm2db.IterationStat.statstype_obj
                ).filter(and_(cm2db.StatsType.category == 'network',
                              cm2db.StatsType.name == 'STRING',
                              cm2db.IterationStat.iteration == iteration)).one().score
        except:
            string_dens = None

        if string_dens is None:
            string_dens = 1.0
        try:
            meme_pval = session.query(cm2db.IterationStat).join(
                cm2db.IterationStat.statstype_obj
                ).filter(and_(cm2db.StatsType.category == 'seqtype',
                              cm2db.StatsType.name == 'upstream',
                              cm2db.IterationStat.iteration == iteration)).one().score
        except:
            # no value found, set meme_pval to 1.0
            meme_pval = 1.0

        last_meme_iteration = get_last_meme_iteration(outdir)
        meme_out = meme_to_str(outdir, last_meme_iteration, cluster)

        try:
            resid = session.query(cm2db.ClusterStat).filter(and_(cm2db.ClusterStat.cluster == cluster,
                                                                 cm2db.ClusterStat.iteration == iteration)).one().residual
        except:
            resid = None
        if resid is None:
            resid = 1.0

        info_line = '"%s"\t%f\t%d\t%f\t"%s"\t%f\t"%s"\n' % (cols_out, string_dens, cluster, meme_pval, meme_out, resid, rows_out)
        if as_binary:
            outfile.write(info_line.encode('utf-8'))
        else:
            outfile.write(info_line)
