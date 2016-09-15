import json
import os
from cmonkey.tools.util import read_ratios
import string

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

def export_to_gaggle_microformats(conn, result_dir, output_dir):
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from row_members')
    iteration = cursor.fetchone()[0]
    cursor.execute('select distinct cluster from row_members where iteration=? order by cluster', [iteration])
    clusters = [row[0] for row in cursor.fetchall()]
    cursor.execute('select organism,species from run_infos')
    orgcode, species = cursor.fetchone()
    templ = string.Template(GAGGLE_TEMPLATE)
    for cluster in clusters:
        # row members
        cursor.execute('select distinct name from row_members rm join row_names rn on rm.order_num=rn.order_num where cluster=? and iteration=?',
                       [cluster, iteration])
        rm_string = '<div class="gaggle-data genes">\n'
        rm_string += '<span class="gaggle-name hidden">Row members</span>\n'
        rm_string += '<span class="gaggle-species">' + species + '</span>\n'
        rm_string += '<div class="gaggle-namelist"><ol>\n'
        for row in cursor.fetchall():
            rm_string += "<li>" + row[0] + "</li>\n"
        rm_string += "</ol></div>"

        # column members
        cursor.execute('select distinct name from column_members cm join column_names cn on cm.order_num=cn.order_num where cluster=? and iteration=?',
                       [cluster, iteration])
        cm_string = '<div class="gaggle-data conditions">\n'
        cm_string += '<span class="gaggle-name hidden">Column members</span>\n'
        cm_string += '<span class="gaggle-species">' + species + '</span>\n'
        cm_string += '<div class="gaggle-namelist"><ol>\n'
        for row in cursor.fetchall():
            cm_string += "<li>" + row[0] + "</li>\n"
        cm_string += "</ol></div>"

        # motifs
        m_string = ''
        cursor.execute('select rowid,motif_num,seqtype from motif_infos where cluster=? and iteration=?',
                       [cluster, iteration])
        motif_infos = [(motif_id, motif_num, seqtype) for motif_id, motif_num, seqtype in cursor.fetchall()]
        for motif_id, motif_num, seqtype in motif_infos:
            m_string += '<div class="gaggle-data motifs">\n'
            m_string += '<span class="gaggle-name">%s motif %d (%s)</span>\n' % (orgcode, motif_num, seqtype)
            m_string += '<span class="gaggle-species">' + species + '</span>\n'
            cursor.execute('select row,a,c,g,t from motif_pssm_rows where motif_info_id=? order by row',
                           [motif_id])
            m_string += '<div class="gaggle-matrix-tsv">\n'
            m_string += "\tPOSITION\tA\tC\tG\tT\n"
            for rownum, a, c, g, t in cursor.fetchall():
                m_string += "%d\t%f\t%f\t%f\t%f\n" % (rownum, a, c, g, t)
            m_string += "</div></div>\n"

        s = templ.substitute(cluster=cluster, row_members=rm_string, column_members=cm_string, motifs=m_string)
        with open(os.path.join(result_dir, 'cluster-%03d.html' % cluster), 'w') as outfile:
            outfile.write(s)


def cluster_expressions_to_json_file(conn, result_dir, output_dir):
    ratios = read_ratios(result_dir)
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from row_members')
    iteration = cursor.fetchone()[0]
    cursor.execute('select distinct cluster from row_members where iteration=?', [iteration])
    clusters = [row[0] for row in cursor.fetchall()]
    result = {}
    for cluster in clusters:
        cursor.execute('select name from row_names rn join row_members rm on rn.order_num=rm.order_num where cluster=? and iteration=?',
                       [cluster, iteration])
        genes = [row[0] for row in cursor.fetchall()]
        cursor.execute('select name from column_names cn join column_members cm on cn.order_num=cm.order_num where cluster=? and iteration=?',
                       [cluster, iteration])
        cluster_conds = [row[0] for row in cursor.fetchall()]
        cluster_data = ratios.loc[genes, cluster_conds]
        values = [{'gene': gene, 'condition': cond, 'value': cluster_data.values[rindex, cindex]}
                  for rindex, gene in enumerate(genes) for cindex, cond in enumerate(cluster_conds)]
        result[str(cluster)] = values

    buffer = json.dumps(result)
    with open(os.path.join(output_dir, 'cluster_expressions.json'), 'w') as out:
        out.write(buffer)


def export_motif_evalues_tsv(conn, result_dir, output_dir):
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    iteration = cursor.fetchone()[0]
    query = """select m1.cluster,m1.seqtype,m1.evalue as evalue1,m2.evalue as evalue2 from (select distinct cluster,seqtype,evalue from motif_infos where iteration=? and motif_num=1) as m1 left outer join (select distinct cluster,seqtype,evalue from motif_infos where iteration=? and motif_num=2) as m2 on m1.cluster=m2.cluster and m1.seqtype=m2.seqtype"""
    cursor.execute(query, [iteration, iteration])
    with open(os.path.join(output_dir, 'motif_evalues.tsv'), 'w') as outfile:
        for row in cursor.fetchall():
            outfile.write("\t".join(map(str, row)) + '\n')

