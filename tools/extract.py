#!/usr/bin/env python
import sqlite3
import sys

def create_database(conn):
    conn.execute('''create table run_infos (start_time timestamp,
                    finish_time timestamp,
                    num_iterations int, last_iteration int,
                    organism text, species text, num_rows int,
                    num_columns int, num_clusters int)''')

    conn.execute('''create table row_names (order_num int, name text)''')
    conn.execute('''create table column_names (order_num int, name text)''')
    conn.execute('''create table row_members (iteration int, cluster int,
                    order_num int)''')
    conn.execute('''create table column_members (iteration int, cluster int,
                    order_num int)''')

    conn.execute('''create table motif_infos (iteration int, cluster int,
                    seqtype text, motif_num int, evalue decimal)''')
    conn.execute('''create table motif_pssm_rows (motif_info_id int,
                    iteration int, row int, a decimal, c decimal, g decimal,
                    t decimal)''')
    conn.execute('''create table meme_motif_sites (motif_info_id int,
                    seq_name text,
                    reverse boolean, start int, pvalue decimal,
                    flank_left text, seq text, flank_right text)''')
    conn.execute('''create table motif_annotations (motif_info_id int,
                    iteration int, gene_num int,
                    position int, reverse boolean, pvalue decimal)''')

        conn.execute('''create table cluster_residuals (iteration int,
                        cluster int, residual decimal)''')

def dbconn(filename):
    return sqlite3.connect(filename, 15, isolation_level='DEFERRED')

def copy_motifs(src, dest, iter):
    cur = src.cursor()
    cur.execute("select rowid,cluster,seqtype,motif_num,evalue from motif_infos where iteration = ?", [iter])

    with dest:
        dcur = dest.cursor()
        scur = src.cursor()
        for motid, cluster, seqtype, motif_num, evalue in cur.fetchall():
            # copy motif info
            dcur.execute("insert into motif_infos (iteration,cluster,seqtype,motif_num,evalue) values (?,?,?,?,?)", [iter, cluster, seqtype, motif_num, evalue])
            motif_info_id = dcur.lastrowid

            # copy pssm
            scur.execute("select row,a,c,g,t from motif_pssm_rows where iteration=? and motif_info_id=?", [iter,motid])
            for row, a, c, g, t in scur.fetchall():
                dcur.execute("insert into motif_pssm_rows (motif_info_id,iteration,row,a,c,g,t) values (?,?,?,?,?,?,?)", [motif_info_id, iter, row, a, c, g, t])

            # copy motif sites
            scur.execute("select seq_name,reverse,start,pvalue,flank_left,seq,flank_right from meme_motif_sites where motif_info_id=?", [motid])
            for seqname,rev,start,pval,flankleft,seq,flankright in scur.fetchall():
                dcur.execute("insert into meme_motif_sites (motif_info_id,seq_name,reverse,start,pvalue,flank_left,seq,flank_right) values (?,?,?,?,?,?,?,?)", [motif_info_id,seqname,rev,start,pval,flankleft,seq,flankright])
            # copy motif annotations
            scur.execute("select gene_num,position,reverse,pvalue from motif_annotations where motif_info_id=? and iteration=?", [motid, iter])
            for genenum,pos,rev,pval in scur.fetchall():
                dcur.execute("insert into motif_annotations (motif_info_id,iteration,gene_num,position,reverse,pvalue) values (?,?,?,?,?,?)", [motif_info_id, iter, genenum, pos, rev, pval])

        dcur.close()
        scur.close()
    cur.close()

def copy_rowcols(src, dest, iter):
    scur = src.cursor()
    # rows
    scur.execute("select order_num, name from row_names order by order_num")
    with dest:
        dcur = dest.cursor()
        for ordernum, name in scur.fetchall():
            dcur.execute("insert into row_names (order_num,name) values (?,?)", [ordernum,name])
        dcur.close()

    scur.execute("select cluster,order_num from row_members where iteration=?", [iter])
    with dest:
        dcur = dest.cursor()
        for cluster, ordernum in scur.fetchall():
            dcur.execute("insert into row_members (iteration,cluster,order_num) values (?,?,?)",
                         [iter,ordernum,name])
        dcur.close()

    # columns
    scur.execute("select order_num, name from column_names order by order_num")
    with dest:
        dcur = dest.cursor()
        for ordernum, name in scur.fetchall():
            dcur.execute("insert into column_names (order_num,name) values (?,?)",
                         [ordernum,name])
        dcur.close()
    
    scur.execute("select cluster,order_num from column_members where iteration=?", [iter])
    with dest:
        dcur = dest.cursor()
        for cluster, ordernum in scur.fetchall():
            dcur.execute("insert into column_members (iteration,cluster,order_num) values (?,?,?)",
                         [iter,ordernum,name])
        dcur.close()
    scur.close()


def copy_essentials(src, dest):
    cur = src.cursor()
    cur.execute('select start_time, finish_time, num_iterations, last_iteration, organism, species, num_rows, num_columns, num_clusters from run_infos')
    stime, ftime, niter, lastiter, organism, species, nrows, ncols, nclust = cur.fetchone()
    with dest:
        dest.execute("insert into run_infos (start_time,finish_time,num_iterations,last_iteration,organism,species,num_rows,num_columns,num_clusters) values(?,?,?,?,?,?,?,?,?)",
                     [stime,ftime,niter,lastiter,organism,species,nrows,ncols,nclust])
    cur.execute('select max(iteration) from row_members')
    lastiter = cur.fetchone()[0]

    print "copying residuals..."
    cur.execute('select cluster, residual from cluster_residuals where iteration=?', [lastiter])
    with dest:
        for cluster, residual in cur.fetchall():
            dest.execute('insert into cluster_residuals (iteration, cluster, residual) values (?,?,?)',
                         [lastiter, cluster, residual])
    cur.close()
    copy_rowcols(src, dest, lastiter)
    copy_motifs(src, dest, lastiter)

    
if __name__ == '__main__':
    src = dbconn(sys.argv[1])
    dest = dbconn(sys.argv[2])
    create_database(dest)
    copy_essentials(src, dest)
    src.close()
    dest.close()
