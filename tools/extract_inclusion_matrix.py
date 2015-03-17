"""extract_includsion_matrix.py - the extract included conditions from an 
output database.  Will create a matrix with cluster numbers as rows and 
experiment names as columns.  Each element will have the value 
'UP', 'DOWN', 'ZERO', or 'EXCLUDED'

This file is part of cMonkey Python. Please see README and LICENSE for
more information and licensing details.
"""
import sys
import sqlite3
import pdb
sys.path.append('../cmonkey/')
import datamatrix as dm
import util

if __name__ == '__main__':
  if len(sys.argv) <= 2:
    print "usage: python %s <out_directory> <outfile.tsv>"
    print "\t example: python %s ../out ../out/inclusionMatrix.tsv"
  else:
    db_file = sys.argv[1]+'/cmonkey_run.db'
    ratio_file = sys.argv[1]+'/ratios.tsv.gz'
    out_file = sys.argv[2]
    
    con = sqlite3.connect(db_file, 15, isolation_level='DEFERRED')
    cur = con.cursor()
    
    cur.execute("select num_clusters from run_infos")
    [num_clusters] = cur.fetchone()
                
    cur.execute("select max(iteration) from cluster_stats")
    [last_iteration] = cur.fetchone()
                
    cur.execute("select name from row_names order by order_num")
    row_names = cur.fetchall()
    row_names = [i[0] for i in row_names]
    
    cur.execute("select name from column_names order by order_num")
    col_names = cur.fetchall()
    col_names = [i[0] for i in col_names]
    
    matrix_factory = dm.DataMatrixFactory([dm.center_scale_filter])
    infile = util.read_dfile(ratio_file, has_header=True, quote='\"')
    ratios = matrix_factory.create_from(infile, False)
    
    #Print out the header line
    f = open(out_file, 'w')
    headerline = 'Cluster'
    for col in col_names:
        headerline = headerline + '\t' + col
    
    f.write(headerline+'\n')
    for cn in range(0,num_clusters):
        #Find the experiments included in this cluster
        cur.execute("select order_num from column_members where iteration = ? and cluster = ?", (last_iteration, cn+1))
        col_members = cur.fetchall()
        col_members = [i[0] for i in col_members]
        col_members.sort()
        cur_names = [col_names[i] for i in col_members]
        
        cur.execute("select order_num from row_members where iteration = ? and cluster = ?", (last_iteration, cn+1))
        row_members = cur.fetchall()
        row_members = [i[0] for i in row_members]
        cur_genes = [row_names[i] for i in row_members]
        
        #Logic: if included, find out if it's up down or zero in ratios.  Otherwise, set it to excluded.
        cur_row = "Cluster%04d" % (cn+1)
        print cur_row
        #pdb.set_trace()
        for col_name in col_names:
            if col_name in cur_names:
                cur_val = ratios.values[row_members,ratios.column_names.index(col_name)].mean()
                if cur_val == 0:
                    cur_row = cur_row + '\tZERO' #Unlikely
                elif cur_val > 0:
                    cur_row = cur_row + '\tUP' 
                else: #cur_val <0
                    cur_row = cur_row + '\tDOWN' 
            else:
                cur_row = cur_row + '\tEXCLUDE'
        f.write(cur_row+'\n')
    f.close()
    cur.close()