import sys
import os
import sqlite3
import numpy as np

class fastaai_db_acc_first_updater:
	def __init__(self, dbpath):
		self.path = dbpath
		self.conn = None
		self.curs = None
		
		self.accessions = None
		
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
			
	def get_accs(self):
		self.accessions = []
		sql = "SELECT scp_acc FROM scp_data"
		for r in self.curs.execute(sql).fetchall():
			self.accessions.append(r[0])
			
	def load_and_convert_table(self, acc):
		print("Working on", acc)
		sql = 'SELECT * FROM "{acc}_genomes"'.format(acc=acc)
		tetra_first_rep = {}
		for r in self.curs.execute(sql).fetchall():
			genome_id = r[0]
			tetras = np.frombuffer(r[1], dtype = np.int32)
			for t in tetras:
				if t not in tetra_first_rep:
					tetra_first_rep[t] = []
				tetra_first_rep[t].append(genome_id)
			
		to_add = []
		for t in sorted(tetra_first_rep):
			addition = np.array(tetra_first_rep[t], dtype = np.int32)
			tetra_first_rep[t] = None
			addition = np.sort(addition)
			next_row = (int(t), addition.tobytes(),)
			#next_row = (int(t), addition,)
			to_add.append(next_row)
		
		print("Adding", len(to_add), "tetramer rows")
		sql = 'INSERT OR REPLACE INTO "{acc}_tetras" VALUES (?, ?)'.format(acc=acc)
		self.curs.executemany(sql, to_add)
	
import argparse
def cleanup_options():
	parser = argparse.ArgumentParser(description='Finalize a FastAAI v2 database. Must be run on each database involved in a query before proceeding.')
	
	parser.add_argument('-db', '--database', dest='database', default = None,
						help='An initialized FastAAI v2 database')	
		
	args = parser.parse_known_args()
	return parser, vars(args[0])
	
	
def finalize_db():
	parser, args = cleanup_options()
	
	if len(sys.argv) < 2:
		parser.print_help()
	else:
		db = args['database']
		if db is not None:
		
			if not os.path.exists(db):
				print("A FastAAI v2 database must already exist at path", db, "or this step cannot proceed.")
				print("Please check that a database already exists or make one with fastaai2 init")
			else:
				mn = fastaai_db_acc_first_updater(db)
				mn.open()
				mn.get_accs()
				for acc in mn.accessions:
					mn.load_and_convert_table(acc)

				mn.conn.commit()
				mn.close()
		else:
			print("Specify a FastAAI V2 database with --database")
			parser.print_help()
		