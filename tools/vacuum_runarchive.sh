#!/bin/bash

# This is a shell script to cleanup a cmonkey-python archive
# It does so by unpacking it, removing unnecessary pkl/runlog files and
# data from non-final iterations out of the database, it then archives
# the result into a new archive
TMP=tmp

echo "cmonkey result cleanup script"
echo "processing file $1..."
FILENAME=`basename $1`
OUTDIR=`echo $FILENAME | cut -d"." -f 1`
if [ ! -d $TMP ]; then
  echo "$TMP does not yet exist, creating"
  mkdir $TMP
fi

echo "unpacking $1..."
tar xfz $1 -C tmp

# 1. check if everything is ok
echo "checking validity of run..."

if [ -f $TMP/$OUTDIR/cmonkey_run.db ]; then
  FINISH=`echo "select finish_time from run_infos;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db`

  if [ -n "$FINISH" ]; then
    echo "Run was finished at: $FINISH"
    echo "removing pkl and runlog files..."
    rm -f $TMP/$OUTDIR/*.pkl
    rm -f $TMP/$OUTDIR/*.runlog
    LAST_ITER=`echo "select max(iteration) from row_members;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db`
    echo "last iteration is: $LAST_ITER"

    echo "remove unused motif information..."
    echo "delete from meme_motif_sites where motif_info_id in (select rowid from motif_infos where iteration < $LAST_ITER);" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "delete from motif_annotations where motif_info_id in (select rowid from motif_infos where iteration < $LAST_ITER);" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "delete from motif_pssm_rows where motif_info_id in (select rowid from motif_infos where iteration < $LAST_ITER);" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "delete from motif_infos where iteration < $LAST_ITER;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "remove unused row/column memberships..."
    echo "delete from row_members where iteration < $LAST_ITER;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "delete from column_members where iteration < $LAST_ITER;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db

    echo "delete statistics..."
    echo "delete from cluster_residuals where iteration < $LAST_ITER;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "delete from motif_pvalues where iteration < $LAST_ITER;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "drop table if exists cluster_stats;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "drop table if exists iteration_stats;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "drop table if exists motif_stats;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db
    echo "drop table if exists network_stats;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db

    echo "vacuuming the database..."
    echo "vacuum;" | sqlite3 $TMP/$OUTDIR/cmonkey_run.db

    echo "repackaging..."
    tar cfz $TMP/$OUTDIR.tar.gz $TMP/$OUTDIR -C $TMP
    echo "Done !"
  else
    echo "ERROR: This run was not finished."
  fi
else
  echo "ERROR: uncompress failed."
fi
