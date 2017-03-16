"""plot_motifs.py - prints the cmonkey2 motifs as PNG images"""
import os
import subprocess
from sqlalchemy import func
import cmonkey.database as cm2db


def make_weblogo(session, motif_id, inpath, outpath):
    """call weblogo on
    weblogo -a "ACGT" -f test.txt -F png -o test.png"""
    try:
        with open(inpath, 'w') as out:
            out.write("\t".join(["PO", "A", "C", "G", "T"]) + "\n")
            for r in session.query(cm2db.MotifPSSMRow).filter(
                    cm2db.MotifPSSMRow.motif_info_id == motif_id).order_by(cm2db.MotifPSSMRow.row):
                row = [r.row, r.a, r.c, r.g, r.t]
                out.write("\t".join(map(str, row)) + "\n")
        subprocess.call(["weblogo", "-a", "ACGT", "-c", "classic", "-f", inpath, "-F", "png", "-o", outpath])
    finally:
        if os.path.exists(inpath):
            os.remove(inpath)


def generate_plots(session, output_dir):
    iteration = session.query(func.max(cm2db.MotifInfo.iteration))
    for m in session.query(cm2db.MotifInfo).filter(cm2db.MotifInfo.iteration == iteration):
        infile = os.path.join(output_dir, "motif_%d.txt" % m.rowid)
        outfile = os.path.join(output_dir, "motif_%d_%d.png" % (m.cluster, m.motif_num))
        make_weblogo(session, m.rowid, infile, outfile)
