import sys
import os
import numpy as np
import sqlite3 

try:
	from hmm_jsonifier import hmm_jsonifier
except:
	from .hmm_jsonifier import hmm_jsonifier

class fastaai_metadata_loader:
	def __init__(self, path):
		self.path = path
		self.conn = None
		self.curs = None
		

		#scp metadata
		self.scp_name_to_id = None
		self.scp_id_to_name = None
		self.scp_id_to_model = None
		self.next_scp_id = 0
		
		#Genome-level data
		self.gix = None
		self.reverse_gix = None
		self.next_genid = 0
		self.protein_counts = None
		self.genome_tax = None
		self.genlens = None
		self.prot_counts = None
		self.genome_classes = None
		self.genome_class_labels = None
		
		self.gen_tax = None
		
		self.gak = None
		self.gas = None
		
		self.recovered_accessions = None
	
	
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
		
		
	def load_scps(self):
		self.scp_name_to_id = {}
		self.scp_id_to_name = {}
		self.scp_id_to_model = {}
		self.next_scp_id = 0
		
		sql = "SELECT * FROM scp_data"
		for result in self.curs.execute(sql).fetchall():
			name = result[1]
			#accession = result[1]
			scpid = result[3]
			#organism_type = result[2]
			scp_model = result[4]
		
			self.scp_name_to_id[name] = scpid
			self.scp_id_to_name[scpid] = name
			self.scp_id_to_model[scpid] = scp_model
			
		self.next_scp_id = len(self.scp_name_to_id)
		
	def load_gix(self):
		self.gix = {}
		self.reverse_gix = {}
		self.protein_counts = {}
		self.genome_tax = {}
		self.genlens = {}
		self.prot_counts = {}
		self.genome_classes = {}
		self.genome_class_labels = {}
		self.gen_tax = {}
		sql = "SELECT * FROM genome_metadata"
		for result in self.curs.execute(sql).fetchall():
			name = result[0]
			genid = result[1]
			genlen = result[2]
			protein_count = result[3]
			scp_count = result[4]
			genome_class = result[5]
			genome_class_label = result[6]
			tax = result[7:]
			self.gix[name] = genid
			self.reverse_gix[genid] = name
			self.protein_counts[genid] = protein_count
			self.genome_tax[genid] = tax
			self.genlens[genid] = genlen
			self.prot_counts[genid] = scp_count
			self.genome_classes[genid] = genome_class
			self.genome_class_labels[genid] = genome_class_label
			
			self.gen_tax[genid] = tax
			
		self.next_genid = len(self.gix)
		
			
	def load_gak(self):
		self.gak = {}
		self.gas = {}
		self.recovered_accessions = []
		
		sql = "SELECT * FROM genome_accession_counts"
		for result in self.curs.execute(sql).fetchall():
			genid = result[0]
			scpid = result[1]
			protein_length = result[2]
			hmm_score = result[3]
			tetra_count = result[4]
			
			self.recovered_accessions.append(scpid)
			
			if genid not in self.gak:
				self.gak[genid] = {}
				self.gas[genid] = {}
				
			self.gak[genid][scpid] = tetra_count
			self.gas[genid][scpid] = hmm_score
	
		recovered_acc_names = []
		self.recovered_accessions = set(self.recovered_accessions)
		for scpid in self.recovered_accessions:
			scpname = self.scp_id_to_name[scpid]
			recovered_acc_names.append(scpname)
			
		self.recovered_accessions = set(recovered_acc_names)

		
	def load_meta(self, skip_gak = False):
		self.open()
		self.load_gix()
		self.load_scps()
		
		if not skip_gak:
			self.load_gak()
			
		self.close()
		
		