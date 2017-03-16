import json
import os
import string
from sqlalchemy import func, and_

from cmonkey.tools.util import read_ratios
import cmonkey.database as cm2db

GAGGLE_TEMPLATE = """
<!doctype html>
<html>
  <head><title>Cluster $cluster Microformats</title></head>
  <body>
    $row_members
    $column_members
    $motifs
  </body>
</html>
"""

def export_to_gaggle_microformats(session, result_dir, output_dir):

    iteration = session.query(func.max(cm2db.RowMember.iteration))
    clusters = [r[0] for r in session.query(cm2db.RowMember.cluster).distinct().filter(
        cm2db.RowMember.iteration == iteration)]
    runinfo = session.query(cm2db.RunInfo).one()
    species = runinfo.species

    templ = string.Template(GAGGLE_TEMPLATE)
    for cluster in clusters:
        # row members
        rm_string = '<div class="gaggle-data genes">\n'
        rm_string += '<span class="gaggle-name hidden">Row members</span>\n'
        rm_string += '<span class="gaggle-species">' + species + '</span>\n'
        rm_string += '<div class="gaggle-namelist"><ol>\n'

        row_names = [r.row_name.name for r in session.query(cm2db.RowMember).filter(
            and_(cm2db.RowMember.cluster == cluster, cm2db.RowMember.iteration == iteration))]
        for row_name in row_names:
            rm_string += "<li>" + row_name + "</li>\n"
        rm_string += "</ol></div>"

        # column members
        cm_string = '<div class="gaggle-data conditions">\n'
        cm_string += '<span class="gaggle-name hidden">Column members</span>\n'
        cm_string += '<span class="gaggle-species">' + species + '</span>\n'
        cm_string += '<div class="gaggle-namelist"><ol>\n'

        cluster_conds = [c.column_name.name for c in session.query(cm2db.ColumnMember).filter(
            and_(cm2db.ColumnMember.cluster == cluster, cm2db.ColumnMember.iteration == iteration))]
        for col_name in cluster_conds:
            cm_string += "<li>" + col_name + "</li>\n"
        cm_string += "</ol></div>"

        # motifs
        m_string = ''

        for m in session.query(cm2db.MotifInfo).filter(
                and_(cm2db.MotifInfo.cluster == cluster, cm2db.MotifInfo.iteration == iteration)):
            m_string += '<div class="gaggle-data motifs">\n'
            m_string += '<span class="gaggle-name">%s motif %d (%s)</span>\n' % (runinfo.organism, m.motif_num, m.seqtype)
            m_string += '<span class="gaggle-species">' + species + '</span>\n'
            m_string += '<div class="gaggle-matrix-tsv">\n'
            m_string += "\tPOSITION\tA\tC\tG\tT\n"

            for r in m.pssm_rows:
                m_string += "%d\t%f\t%f\t%f\t%f\n" % (r.row, r.a, r.c, r.g, r.t)
            m_string += "</div></div>\n"

        s = templ.substitute(cluster=cluster, row_members=rm_string, column_members=cm_string, motifs=m_string)
        with open(os.path.join(output_dir, 'cluster-%03d.html' % cluster), 'w') as outfile:
            outfile.write(s)


def cluster_expressions_to_json_file(session, result_dir, output_dir):
    ratios = read_ratios(result_dir)

    iteration = session.query(func.max(cm2db.RowMember.iteration))
    clusters = [r[0] for r in session.query(cm2db.RowMember.cluster).distinct().filter(
        cm2db.RowMember.iteration == iteration)]

    result = {}
    for cluster in clusters:
        genes = [r.row_name.name for r in session.query(cm2db.RowMember).filter(
            and_(cm2db.RowMember.cluster == cluster, cm2db.RowMember.iteration == iteration))]
        cluster_conds = [c.column_name.name for c in session.query(cm2db.ColumnMember).filter(
            and_(cm2db.ColumnMember.cluster == cluster, cm2db.ColumnMember.iteration == iteration))]
        cluster_data = ratios.loc[genes, cluster_conds]
        values = [{'gene': gene, 'condition': cond, 'value': cluster_data.values[rindex, cindex]}
                  for rindex, gene in enumerate(genes) for cindex, cond in enumerate(cluster_conds)]
        result[str(cluster)] = values

    buffer = json.dumps(result)
    with open(os.path.join(output_dir, 'cluster_expressions.json'), 'w') as out:
        out.write(buffer)


def export_motif_evalues_tsv(session, result_dir, output_dir):
    iteration = session.query(func.max(cm2db.MotifInfo.iteration))
    clusters = [r[0] for r in session.query(cm2db.MotifInfo.cluster).distinct().filter(
        cm2db.MotifInfo.iteration == iteration)]
    with open(os.path.join(output_dir, 'motif_evalues.tsv'), 'w') as outfile:
        for cluster in clusters:
            motifs = list(session.query(cm2db.MotifInfo).filter(
                and_(cm2db.MotifInfo.iteration == iteration, cm2db.MotifInfo.cluster == cluster)))
            if len(motifs) == 2:
                row = [str(cluster), motifs[0].seqtype, str(motifs[0].evalue), str(motifs[1].evalue)]
                outline = '\t'.join(row) + '\n'
                outfile.write(outline)

