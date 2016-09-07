"""plot_motif_positions.py - generate motif annotation plots"""
import os
import svgwrite


def generate_plots(conn, output_dir):
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    max_iteration = cursor.fetchone()[0]
    cursor.execute('select distinct cluster from motif_infos where iteration=?',
                   [max_iteration])
    for row in cursor.fetchall():
        cluster = row[0]
        dwg = svgwrite.Drawing(os.path.join(output_dir, 'motif_pos-%d.svg' % cluster))
        dwg.add(dwg.line((10, 10), (20, 10), stroke='black'))
        dwg.save()
        break
