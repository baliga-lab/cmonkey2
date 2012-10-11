"""sequence_cache.py - sqlite-based cache for genes sequences
This module supports offloading gene sequences into a sqlite database.
cMonkey rarely requires to access all sequences, but storing sequence
information in memory potentially can take up a huge amount of space.
The database cached to support sequence types with both search and
scan distances.

The database is optimized to support the following queries:

- by sequence type and gene name
- by sequence type
"""
import sqlite3
import os.path

class SequenceCache:
    """schema
    sequence_type(id, name, start, stop)
    sequence(id, sequence_type_id, sequence_name, sequence)
    """
    def __init__(self, path, create=True):
        """creates the cache from the specified path. If the database file
        exists, it will be deleted
        """
        if create:
            if os.path.exists(path):
                os.remove(path)

        self.conn = sqlite3.connect(path)

        if create:
            c = self.conn.cursor()
            c.execute('''create table sequence_types
                         (name text, start int, stop int)''')
            c.execute('''create table sequences
                         (sequence_type_id int, name text, contig text,
                          start int, end int, is_reverse int, sequence text)''')
            c.execute('''create index seq_nt_index
                         on sequences (sequence_type_id, name)''')
            c.execute('''create index seq_type_index
                         on sequences (sequence_type_id)''')
            self.conn.commit()
            c.close()

    def add_sequence_type(self, name, start, stop):
        c = self.conn.cursor()
        c.execute('''insert into sequence_types (name, start, stop)
                     values (?, ?, ?)''', (name, start, stop))
        self.conn.commit()
        c.close()

    def add_sequence(self, seqtype_name, start, stop, seq_name, sequence,
                     contig, loc_start, loc_end, is_reverse):
        c = self.conn.cursor()
        c.execute('''select ROWID from sequence_types where
                     name = ? and start = ? and stop = ?''', (seqtype_name, start, stop))
        rowid = c.fetchone()[0]
        c.execute('''insert into sequences (sequence_type_id, name, sequence, contig,
                                            start, end, is_reverse)
                     values (?, ?, ?, ?, ?, ?, ?)''', (rowid, seq_name, sequence, contig,
                                                    loc_start, loc_end, is_reverse))
        self.conn.commit()
        c.close()

    def add_sequences(self, seqtype_name, start, stop, data):
        """Adding data all at once, which should be faster when inserting a lot
        of sequences.
        data must be of format (gene, sequence, contig, loc_start, loc_end, reverse)"""
        c = self.conn.cursor()
        c.execute('''select ROWID from sequence_types where
                     name = ? and start = ? and stop = ?''', (seqtype_name, start, stop))
        rowid = c.fetchone()[0]
        records = [([rowid] + list(row)) for row in data]
        c.executemany('''insert into sequences (sequence_type_id, name, sequence, contig,
                                            start, end, is_reverse)
                         values (?, ?, ?, ?, ?, ?, ?)''', records)
        self.conn.commit()
        c.close()
        

    def sequences_for_type(self, seqtype_name, start, stop):
        c = self.conn.cursor()
        c.execute('''select ROWID from sequence_types where
                     name = ? and start = ? and stop = ?''', (seqtype_name, start, stop))
        rowid = c.fetchone()[0]
        c.execute('''select name, sequence from sequences where
                     sequence_type_id = ?''', (rowid,))
        results = c.fetchall()
        c.close()
        return results

    def sequences_for_type_and_name(self, seqtype_name, start, stop, seq_names):
        c = self.conn.cursor()
        c.execute('''select ROWID from sequence_types where
                     name = ? and start = ? and stop = ?''', (seqtype_name, start, stop))
        print "seqtype: ", seqtype_name, " start: ", start, " stop: ", stop
        rowid = c.fetchone()[0]
        params = [rowid]
        params.extend(seq_names)
        query = '''select name, sequence from sequences where
                   sequence_type_id = ? and name in (%s)''' % ','.join('?'*len(seq_names))
        c.execute(query, params)
        results = c.fetchall()
        c.close()
        return [(row[0], row[1]) for row in results]
