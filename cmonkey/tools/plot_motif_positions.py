"""plot_motif_positions.py - generate motif annotation plots"""
import os
import svgwrite
from collections import defaultdict
from sqlalchemy import func, and_
import cmonkey.database as cm2db

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


def draw_annotations(session, output_dir, motif_lengths, iteration, cluster):

    # TODO: for multiple seqtypes, we need to outfactor the query to generate a graph
    # for each seqtype
    motif_infos = session.query(cm2db.MotifInfo).filter(
        and_(cm2db.MotifInfo.iteration == iteration, cm2db.MotifInfo.cluster == cluster))

    st_annots = defaultdict(list)  # group by seqtype
    for m in motif_infos:
        for a in m.annotations:
            annot = (m.rowid, m.seqtype, m.motif_num, a.gene.name, a.position, a.reverse, a.pvalue)
            st_annots[m.seqtype].append(annot)

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


def generate_plots(session, output_dir):
    iteration = session.query(func.max(cm2db.RowMember.iteration))
    motif_lengths = {m_id: nrows for m_id, nrows in session.query(
        cm2db.MotifPSSMRow.motif_info_id, func.count(cm2db.MotifPSSMRow.row)).filter(
            cm2db.MotifPSSMRow.iteration == iteration).group_by(cm2db.MotifPSSMRow.motif_info_id)}
    for row in session.query(cm2db.MotifInfo.cluster).distinct().filter(cm2db.MotifInfo.iteration == iteration):
        draw_annotations(session, output_dir, motif_lengths, iteration, row[0])
