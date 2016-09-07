"""plot_expressions.py - make cluster gene expression plots"""
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from cmonkey.tools.util import read_ratios


def normalize_js(value):
    if math.isnan(value) or math.isinf(value):
        return 0.0
    else:
        return value


def generate_plots(conn, result_dir, output_dir):
    ratios = read_ratios(result_dir)

    cursor = conn.cursor()
    cursor.execute('select max(iteration) from row_members')
    iteration = cursor.fetchone()[0]

    cursor.execute('select distinct cluster from row_members where iteration=?', [iteration])
    clusters = [row[0] for row in cursor.fetchall()]
    for cluster in clusters:
        plt.clf()
        plt.cla()
        plt.figure(figsize=(6,3))
        cursor.execute('select distinct name from row_members rm join row_names rn on rm.order_num=rn.order_num where cluster=? and iteration=?', [cluster, iteration])
        genes = [row[0] for row in cursor.fetchall()]

        cursor.execute('select distinct name from column_members cm join column_names cn on cm.order_num=cn.order_num where cluster=? and iteration=?', [cluster, iteration])
        cluster_conds = [row[0] for row in cursor.fetchall()]

        cursor.execute('select distinct name from column_names')
        all_conds = [row[0] for row in cursor.fetchall()]
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
