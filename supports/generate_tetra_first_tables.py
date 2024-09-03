import sys
import os
import sqlite3
import numpy as np

from metadata_loader import fastaai_metadata_loader

import tempfile

import multiprocessing

class partial_fastaai_db:
	def __init__(self, dbpath):
		self.path = dbpath
		self.conn = None
		self.curs = None
		
		self.metadata = None
		self.next_genome_id = 0
		self.next_scp_id = 0
		
		self.accession_index = None
		self.genome_index = None
		self.reverse_acc = None

		
	def open(self):
		if self.conn is None:
			self.conn = sqlite3.connect(self.path)
			self.curs = self.conn.cursor()
	
	def close(self):
		if self.conn is not None:
			self.curs.close()
			self.conn.close()
			self.conn = None
			self.curs = None
			
	def get_meta(self):
		self.metadata = fastaai_metadata_loader(self.path)
		self.metadata.load_meta(skip_gak = False)
		self.next_genome_id = self.metadata.next_genid
		self.next_scp_id = self.metadata.next_scp_id
		
		self.accession_index = self.metadata.scp_name_to_id
		self.genome_index = self.metadata.gix
		self.reverse_acc = self.metadata.scp_id_to_name
		
	def prep_tetra_first_table(self, scp_label):
		sql = 'DROP TABLE IF EXISTS "{scp}_tetras"'.format(scp = scp_label)
		
		tetra_first = ', '.join([
				"tetramer INTEGER",
				"genomes BLOB"
				])
				
		sql = 'CREATE TABLE IF NOT EXISTS "{scp}_tetras" ({fields})'.format(scp = scp_label, fields = tetra_first)
		self.curs.execute(sql)

		sql = 'DROP INDEX IF EXISTS "{scp}_tetra_first"'

		sql = 'CREATE INDEX IF NOT EXISTS "{scp}_tetra_first" ON "{scp}_tetras" (tetramer)'.format(scp = scp_label)
		self.curs.execute(sql)
		
	
def converter_init(refdb_path):
	global refdb
	refdb = partial_fastaai_db(refdb_path)
	
	return refdb
	
def load_and_convert_one_scp(scp):
	#print("Working on", scp)
	refdb.open()
	sql = 'SELECT * FROM "{scp}_genomes"'.format(scp = scp)
	results = refdb.curs.execute(sql).fetchall()
	refdb.close()

	#print(scp, "retrieved", len(results), "genomes")

	reformatted = {}
	for r in results:
		genome, tetras = r[0], np.frombuffer(r[1], dtype = np.int32)
		for t in tetras:
			if t not in reformatted:
				reformatted[t] = []
			reformatted[t].append(genome)
		
	results = None
		
	insertable = []
	for t in reformatted:
		arr = np.array(reformatted[t], dtype = np.int32)
		arr = np.sort(arr)
		next_row = (t, arr,)
		insertable.append(next_row)
		reformatted[t] = None
		
	reformatted = None	

	#print(scp, "converted to", len(insertable), "tetramers")

	tmp = tempfile.TemporaryDirectory()
	tmppath = tmp.name
	tmppath = os.path.normpath(tmppath+"_"+scp+"_.sqlite3.db")
	tempdb = tmppath
	
	mn = partial_fastaai_db(tempdb)
	mn.open()
	mn.prep_tetra_first_table(scp)
	mn.curs.executemany('INSERT INTO "{scp}_tetras" VALUES (?, ?)'.format(scp = scp), insertable)
	mn.conn.commit()
	insertable = None
	mn.close()

	return scp, tmppath

	
f = "comprehensive_proks.db"
	
mn = partial_fastaai_db(f)
mn.open()
mn.get_meta()
mn.close()

scps = list(mn.accession_index.keys())
	

thds = int(sys.argv[1])


results = []
pool = multiprocessing.Pool(thds, initializer = converter_init, initargs = (f,))

for scp, dbpath in pool.imap_unordered(load_and_convert_one_scp, scps):
	results.append((scp, dbpath,))
	print(scp, "reordering generated!")
	
pool.close()
pool.join()


mn.open()
for pair in results:
	scp = pair[0]
	dbpath = pair[1]
	
	mn.prep_tetra_first_table(scp)
	
	print(scp, "preparing for final database insert...")
	sub = partial_fastaai_db(dbpath)
	sub.open()
	results = sub.curs.execute('SELECT * FROM "{scp}_tetras"'.format(scp = scp)).fetchall()
	sub.close()
	
	sql = 'INSERT INTO "{scp}_tetras" VALUES (?, ?)'.format(scp = scp)
	mn.curs.executemany(sql, results)
	mn.conn.commit()
	
	print(scp, "inserted into final database!")
	os.remove(dbpath)
	

mn.close()