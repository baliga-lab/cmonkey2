import json
import os
from cmonkey.tools.util import read_ratios


def to_json_file(conn, result_dir, output_dir):
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

