"""plot_motif_positions.py - generate motif annotation plots"""
import os
import svgwrite
from collections import defaultdict

TEXT_STYLE = "font-size:%ipx; font-family:%s" % (10, "sans-serif")
BOX_WIDTH = 100.0
BOX_HEIGHT = 20.0
BOX_VSPACE = 25.0
MARGIN_RIGHT = 20.0
MARGIN_TOP   = 20.0
BASE_WIDTH = 500.0
SCALE_WIDTH = 200.0
TICK_LEN = 7.0
SCALE_OFFSET = 30
BASE_X_OFFSET = 140.0
MATCH_HEIGHT = 10.0
PVALUE_CUTOFF = 0.6


def draw_scale(dwg, base_x, height):
    tick_labels = ['-200', '-100', '-1']

    x = base_x - SCALE_WIDTH
    y = height - SCALE_OFFSET
    right = x + SCALE_WIDTH
    tick_top = y - TICK_LEN
    middle = x + SCALE_WIDTH / 2.0

    def draw_tick(x, label):
        dwg.add(dwg.line((x, tick_top), (x, y), stroke='black'))
        text = dwg.text(label, x=[x - 8], y=[y + 15], stroke='black', style=TEXT_STYLE)
        dwg.add(text)

    dwg.add(dwg.line((x, y), (right, y), stroke='black'))
    draw_tick(x, tick_labels[0])
    draw_tick(middle, tick_labels[1])
    draw_tick(right, tick_labels[2])


def draw_annotation(dwg, annotation, base_x, annot_y, motif_lengths):
    name, color, matches = annotation
    box_y = annot_y
    text_offset_y = BOX_HEIGHT * 0.75
    box_x = dwg.attribs['width'] - BOX_WIDTH - MARGIN_RIGHT
    line_x = 40
    line_y = box_y + BOX_HEIGHT / 2.0

    dwg.add(dwg.rect((box_x, box_y), (BOX_WIDTH, BOX_HEIGHT), stroke=color, fill=color))
    dwg.add(dwg.line((line_x, line_y), (box_x, line_y), stroke='gray'))
    text = dwg.text(name, x=[box_x + BOX_WIDTH / 2.0], y=[box_y + text_offset_y], stroke='black', style=TEXT_STYLE)
    dwg.add(text)

    for motif_id, seqtype, motif_num, name, position, reverse, pvalue in matches:
        mlen = motif_lengths[motif_id]
        mx = base_x - position - mlen
        my = line_y if reverse else line_y - MATCH_HEIGHT
        mcolor = '#f00' if motif_num == 1 else '#0f0'
        mopacity = 0.0 if pvalue > PVALUE_CUTOFF else 1.0 - (pvalue * 1.5)
        dwg.add(dwg.rect((mx, my), (mlen, MATCH_HEIGHT), stroke="none", fill=mcolor, fill_opacity=mopacity))


def draw_annotations(conn, output_dir, motif_lengths, iteration, cluster):

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

        width = BASE_WIDTH
        height = MARGIN_TOP * 2 + len(gene_annots) * BOX_VSPACE + SCALE_OFFSET
        dwg = svgwrite.Drawing(os.path.join(output_dir, 'motif_pos_%s-%d.svg' % (seqtype, cluster)), (width, height))

        # border
        dwg.add(dwg.rect((0, 0), (width, height), stroke="blue", fill="white"))
        base_x = width - BASE_X_OFFSET

        draw_scale(dwg, base_x, height)
        dwg.add(dwg.line((base_x, 0), (base_x, height), stroke='gray', stroke_dasharray=[1, 3]))

        annot_y = MARGIN_TOP
        for gene, matches in gene_annots.items():
            annotation = (gene, '#00ff88', matches)
            draw_annotation(dwg, annotation, base_x, annot_y, motif_lengths)
            annot_y += BOX_VSPACE

        dwg.save()



def generate_plots(conn, output_dir):
    cursor = conn.cursor()
    cursor.execute('select max(iteration) from motif_infos')
    max_iteration = cursor.fetchone()[0]
    cursor.execute('select motif_info_id, count(row) from motif_pssm_rows where iteration = ? group by motif_info_id',
                   [max_iteration])
    motif_lengths = {row[0]: row[1] for row in cursor.fetchall()}

    cursor.execute('select distinct cluster from motif_infos where iteration=?',
                   [max_iteration])
    for row in cursor.fetchall():
        draw_annotations(conn, output_dir, motif_lengths, max_iteration, row[0])
