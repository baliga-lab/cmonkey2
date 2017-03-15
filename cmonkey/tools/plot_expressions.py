"""plot_expressions.py - make cluster gene expression plots"""
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from cmonkey.tools.util import read_ratios
import cmonkey.database as cm2db
from sqlalchemy import func, and_


def normalize_js(value):
    if math.isnan(value) or math.isinf(value):
        return 0.0
    else:
        return value


def generate_plots(session, result_dir, output_dir):
    ratios = read_ratios(result_dir)

    iteration = session.query(func.max(cm2db.RowMember.iteration))
    clusters = [r[0] for r in session.query(cm2db.RowMember.cluster).distinct().filter(
        cm2db.RowMember.iteration == iteration)]

    figure = plt.figure(figsize=(6,3))
    for cluster in clusters:
        plt.clf()
        plt.cla()
        genes = [r.row_name.name for r in session.query(cm2db.RowMember).filter(
            and_(cm2db.RowMember.cluster == cluster, cm2db.RowMember.iteration == iteration))]
        cluster_conds = [c.column_name.name for c in session.query(cm2db.ColumnMember).filter(
            and_(cm2db.ColumnMember.cluster == cluster, cm2db.ColumnMember.iteration == iteration))]
        all_conds = [c[0] for c in session.query(cm2db.ColumnName.name).distinct()]
        non_cluster_conds = [cond for cond in all_conds if not cond in set(cluster_conds)]

        cluster_data = ratios.loc[genes, cluster_conds]
        non_cluster_data = ratios.loc[genes, non_cluster_conds]
        min_value = ratios.min()
        max_value = ratios.max()
        for gene in genes:
            values = [normalize_js(val) for val in cluster_data.loc[gene,:].values]
            values += [normalize_js(val) for val in non_cluster_data.loc[gene,:].values]
            plt.plot(values)

        # plot the "in"/"out" separator line
        cut_line = len(cluster_conds)
        plt.plot([cut_line, cut_line], [min_value, max_value], color='red',
                 linestyle='--', linewidth=1)
        plt.savefig(os.path.join(output_dir, "exp-%d" % cluster))
    plt.close(figure)
