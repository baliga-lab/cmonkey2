#!/bin/bash

# This is a shell script to cleanup a cmonkey-python result archive or directory
# If directory:
#   - removing unnecessary pkl/runlog files and data from database
#
# If it is an archive:
#   - unpacking it
#   - cleanup the unarchived directory
#   - re-archives result into a new archive
TMP=tmp
SCRIPTDIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

function cleanup_outdir {
    if [ -f $1/cmonkey_run.db ]; then
        FINISH=`echo "select finish_time from run_infos;" | sqlite3 $1/cmonkey_run.db`

        if [ -n "$FINISH" ]; then
            echo "Run was finished at: $FINISH"
            echo "removing pkl and runlog files..."
            rm -f $1/*.pkl
            rm -f $1/*.runlog
            echo "compress the database..."
            $SCRIPTDIR/extract.py $1/cmonkey_run.db $1/cmonkey_run_extract.db && mv $1/cmonkey_run_extract.db $1/cmonkey_run.db
            return 0;
        else
            echo "ERROR: This run was not finished."
        fi
    else
        echo "ERROR: uncompress failed."
    fi
    return 1;
}

echo "cmonkey result cleanup script"

if [ -d $1 ]; then
    echo "processing output directory $1..."
    cleanup_outdir $1
else
    echo "processing file $1..."
    FILENAME=`basename $1`
    OUTDIR=`echo $FILENAME | cut -d"." -f 1`
    if [ ! -d $TMP ]; then
        echo "$TMP does not yet exist, creating"
        mkdir $TMP
    fi

    echo "unpacking $1..."
    tar xfz $1 -C $TMP
    if cleanup_outdir "$TMP/$OUTDIR"; then
        echo "repackaging $TMP/$OUTDIR..."
        tar cfz $TMP/$OUTDIR.tar.gz $TMP/$OUTDIR -C $TMP
    fi
fi
echo "Done !"

