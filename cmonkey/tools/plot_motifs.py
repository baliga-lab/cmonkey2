"""plot_motifs.py - prints the cmonkey2 motifs as PNG images"""
import os
import subprocess


def make_weblogo(cursor, motif_id, inpath, outpath):
    """call weblogo on
    weblogo -a "ACGT" -f test.txt -F png -o test.png"""
    try:
        with open(inpath, 'w') as out:
            out.write("\t".join(["PO", "A", "C", "G", "T"]) + "\n")
            cursor.execute('select row,a,c,g,t from motif_pssm_rows where motif_info_id=? order by row',
                           [motif_id])
            for row in cursor.fetchall():
                out.write("\t".join(map(str, row)) + "\n")
        subprocess.call(["weblogo", "-a", "ACGT", "-c", "classic", "-f", inpath, "-F", "png", "-o", outpath])
    finally:
        if os.path.exists(inpath):
            os.remove(inpath)


def generate_plots(conn, output_dir):
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    max_iteration = cursor.fetchone()[0]
    cursor.execute('select rowid,cluster,motif_num from motif_infos where iteration=?',
                   [max_iteration])
    for motif_id, cluster, motif_num in cursor.fetchall():
        infile = os.path.join(output_dir, "motif_%d.txt" % motif_id)
        outfile = os.path.join(output_dir, "motif_%d_%d.png" % (cluster, motif_num))
        make_weblogo(cursor, motif_id, infile, outfile)
