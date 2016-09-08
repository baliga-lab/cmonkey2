"""plot_motif_positions.py - generate motif annotation plots"""
import os
import svgwrite
from collections import defaultdict

TEXT_STYLE = "font-size:%ipx; font-family:%s" % (10, "sans-serif")
BOX_WIDTH = 100.0
BOX_HEIGHT = 20.0
MARGIN_RIGHT = 20.0

def draw_scale(dwg, base_x, height):
    tick_labels = ['-200', '-100', '-1']
    tick_len = 7
    scale_width = 200

    x = base_x - scale_width
    y = height - 30
    right = x + scale_width
    tick_top = y - tick_len
    middle = x + scale_width / 2.0
    def draw_tick(x, label):
        dwg.add(dwg.line((x, tick_top), (x, y), stroke='black'))
        text = dwg.text(label, x=[x - 8], y=[y + 15], stroke='black', style=TEXT_STYLE)
        dwg.add(text)

    dwg.add(dwg.line((x, y), (right, y), stroke='black'))
    draw_tick(x, tick_labels[0])
    draw_tick(middle, tick_labels[1])
    draw_tick(right, tick_labels[2])


def draw_annotation(dwg, annotation, annot_y):
    name, color = annotation
    box_y = annot_y
    text_offset_y = BOX_HEIGHT / 2.0
    box_x = dwg.attribs['width'] - BOX_WIDTH - MARGIN_RIGHT
    line_x = 40
    line_y = box_y + BOX_HEIGHT / 2.0

    dwg.add(dwg.rect((box_x, box_y), (BOX_WIDTH, BOX_HEIGHT), stroke=color, fill=color))
    dwg.add(dwg.line((line_x, line_y), (box_x, line_y), stroke='gray'))
    text = dwg.text(name, x=[box_x + BOX_WIDTH / 2.0], y=[box_y + text_offset_y], stroke='black', style=TEXT_STYLE)
    dwg.add(text)


def draw_annotations(conn, dwg, motif_lengths, iteration, cluster):
    annot_y = 20.0

    # TODO: for multiple seqtypes, we need to outfactor the query to generate a graph
    # for each seqtype
    cursor = conn.cursor()
    cursor.execute('select a.motif_info_id, seqtype, motif_num, g.name, position, reverse, pvalue from motif_annotations a join motif_infos i on a.motif_info_id = i.rowid join row_names g on g.order_num = a.gene_num where i.iteration = ? and i.cluster = ?', [iteration, cluster])

    #annotations = [("gene 1", "#00ff88"), ("gene 2", "#ff1245")]
    annotations = []
    st_annots = defaultdict(list)  # group by seqtype
    for motif_id, seqtype, motif_num, name, position, reverse, pvalue in cursor.fetchall():
        annot = (motif_id, seqtype, motif_num, name, position, reverse, pvalue)
        st_annots[seqtype].append(annot)

    st_gene_annots = {}
    for seqtype in st_annots:  # group by gene
        gene_annots = defaultdict(list)
        for annot in st_annots[seqtype]:
            gene_annots[annot[3]].append(annot)  # gene name
        for gene, annots in gene_annots.items():
            annotation = (gene, '#00ff88')
            draw_annotation(dwg, annotation, annot_y)
            annot_y += 25.0


def generate_plots(conn, output_dir):
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    max_iteration = cursor.fetchone()[0]
    cursor.execute('select motif_info_id, count(row) from motif_pssm_rows where iteration = ? group by motif_info_id',
                   [max_iteration])
    motif_lengths = {row[0]: row[1] for row in cursor.fetchall()}

    width = 500.0

    cursor.execute('select distinct cluster from motif_infos where iteration=?',
                   [max_iteration])
    for row in cursor.fetchall():
        cluster = row[0]
        height = 400.0
        dwg = svgwrite.Drawing(os.path.join(output_dir, 'motif_pos-%d.svg' % cluster), (width, height))
        # border
        dwg.add(dwg.rect((0, 0), (width, height), stroke="blue", fill="white"))
        base_x = width - 140.0

        draw_scale(dwg, base_x, height)
        dwg.add(dwg.line((base_x, 0), (base_x, height), stroke='gray', stroke_dasharray=[1, 3]))
        draw_annotations(conn, dwg, motif_lengths, max_iteration, cluster)
        dwg.save()
        break
