import sys
import os
import numpy as np
import sqlite3 
import json


from supports.metadata_loader import fastaai_metadata_loader
from supports.crystal_loader import ravenous_crystal_lizard

class fastaai_consumer:
	def __init__(self, dbpath):
		self.path = dbpath
		self.conn = None
		self.curs = None
		
		self.add_rule = None
		
		self.metadata = None
		self.next_genome_id = 0
		self.next_scp_id = 0
		
		self.accession_index = None
		self.genome_index = None
		self.reverse_acc = None
		
		self.condition = None
		self.rule_rider = None
		self.genome_class_name_to_num = None
		
		
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
		'''
		#metadata contents
		#scp metadata
		self.scp_name_to_id = None
		self.scp_id_to_name = None
		self.scp_id_to_model = None
		self.next_scp_id = 0
		
		#Genome-level data
		self.gix = None
		self.reverse_gix = None
		self.next_genid = 0
		self.gen_classes = None
		self.genome_tax = None
		self.genlens = None
		self.prot_counts = None
		
		self.gak = None
		'''
		self.metadata = fastaai_metadata_loader(self.path)
		self.metadata.load_meta(skip_gak = False)
		self.next_genome_id = self.metadata.next_genid
		self.next_scp_id = self.metadata.next_scp_id
		
		'''
		self.scp_name_to_id
		self.scp_id_to_name
		self.scp_id_to_model
		'''
		self.accession_index = self.metadata.scp_name_to_id
		self.genome_index = self.metadata.gix
		self.reverse_acc = self.metadata.scp_id_to_name
		
		#Genome addition rules
		rule = self.curs.execute("SELECT ruleset_json FROM rule_addition_fields").fetchone()
		try:
			rule = rule[0]
			rule = json.loads(rule)
			
		except:
			print("Rule failed to load")
			rule = None
				
		self.add_rule = None
		self.condition = None
		
		if rule is not None:
			if 'class_assignment' in rule:
				self.condition = rule['class_assignment']['add_condition']
				if self.condition == "majority_vote":
					self.genome_class_name_to_num = {}
					self.add_rule = {}
					for gen_class in rule['class_assignment']['classes']:
						self.genome_class_name_to_num[gen_class] = rule['class_assignment']['classes'][gen_class]['class']
						accs_to_add = []
						for acc in rule['class_assignment']['classes'][gen_class]['members']:
							if acc in self.accession_index:
								accession_id = self.accession_index[acc]
								accs_to_add.append(accession_id)
							else:
								print("accession not found:", acc)
					
						self.add_rule[gen_class] = set(accs_to_add)
					
		self.rule_rider = None
		if 'class_assignment' in rule:
			self.rule_rider = rule['class_assignment']['rider']
		else:
			self.rule_rider = None
	


	def get_majority_vote(self, gen_md, gak_dat, genome_first_rep):
		genome = gen_md[0]

		accessions = []
		for g in gak_dat:
			accession_id = g[1]
			accessions.append(accession_id)
			#acc_name = self.reverse_acc[accession_id]
		accessions = set(accessions)
		
		winning_class = None
		winning_score = 0
		
		score_by_class = {}
		for genome_class in self.add_rule:
			score = len(accessions.intersection(self.add_rule[genome_class])) / len(self.add_rule[genome_class])
			score_by_class[genome_class] = score
			if score > winning_score:
				winning_score = score
				winning_class = genome_class
				
		if self.rule_rider == "exclusive":
			acceptable_accessions = self.add_rule[winning_class]
			cleaned_accessions = {}
			for acc in genome_first_rep:
				acc_id = self.accession_index[acc]
				if acc_id in acceptable_accessions:
					cleaned_accessions[acc] = genome_first_rep[acc]
					genome_first_rep[acc] = None
			
			genome_first_rep = cleaned_accessions
				
		count_of_scps = len(genome_first_rep)
		
		winning_class_id = self.genome_class_name_to_num[winning_class]
		
		gen_md[6] = winning_class
		gen_md[5] = winning_class_id
		gen_md[4] = count_of_scps
		
		if self.rule_rider == "exclusive":
			cleaned_gak = []
			for g in gak_dat:
				accession_id = g[1]
				if accession_id in acceptable_accessions:
					cleaned_gak.append(g)
			gak_dat = cleaned_gak
				
		return gen_md, gak_dat, genome_first_rep
		
	
	def add_crystals(self, crystal_list):
		crystal_eater = ravenous_crystal_lizard()
		crystal_eater.set_indices(self.genome_index, self.accession_index, self.next_genome_id)
		
		gix_additions = []
		gak_additions = []
		genome_first_additions = {}
		
		to_add = 0
		for f in crystal_list:
			#Corresponds to tables
			#genome_metadata, genome_accession_counts, {scp}_genome
			gen_md, gak_dat, genome_first_rep = crystal_eater.load_crystal(f)
			
			if self.condition is not None:
				if self.condition == "majority_vote":
					#Filter results as needed
					gen_md, gak_dat, genome_first_rep = self.get_majority_vote(gen_md, gak_dat, genome_first_rep)

			#Probably need to check for a repeat genome add here.
			gix_additions.append(gen_md)
			gak_additions.extend(gak_dat)
			
			for acc in genome_first_rep:
				if acc not in genome_first_additions:
					genome_first_additions[acc] = []
				genome_first_additions[acc].append(genome_first_rep[acc])
			
			to_add += 1
			
			if to_add % 1000 == 0:
				print("Committing next 1000 genomes. Total added:", to_add)
				self.curs.executemany("INSERT OR REPLACE INTO genome_metadata VALUES ({gm_fmt})".format(gm_fmt = ', '.join(['?']*14)), gix_additions)
				self.curs.executemany("INSERT OR REPLACE INTO genome_accession_counts VALUES ({gak_fmt})".format(gak_fmt = ', '.join(['?']*5)), gak_additions)
				for acc in genome_first_additions:
					self.curs.executemany('INSERT OR REPLACE INTO "{acc}_genomes" VALUES (?, ?)'.format(acc = acc), genome_first_additions[acc])
				
				gix_additions = []
				gak_additions = []
				genome_first_additions = {}
				self.conn.commit()
		
		print("Committing final batch of genomes. Total added:", to_add)
		self.curs.executemany("INSERT OR REPLACE INTO genome_metadata VALUES ({gm_fmt})".format(gm_fmt = ', '.join(['?']*14)), gix_additions)
		self.curs.executemany("INSERT OR REPLACE INTO genome_accession_counts VALUES ({gak_fmt})".format(gak_fmt = ', '.join(['?']*5)), gak_additions)
		for acc in genome_first_additions:
			self.curs.executemany('INSERT OR REPLACE INTO "{acc}_genomes" VALUES (?, ?)'.format(acc = acc), genome_first_additions[acc])
			
		gix_additions = []
		gak_additions = []
		genome_first_additions = {}
		self.conn.commit()

def fastaai_consume():
	parser, args = consumer_options()
	
	if len(sys.argv) < 3:
		parser.print_help()
	else:
		crystal_dir = args['crystals_dir']
		dbpath = args['database']
		
		if crystal_dir is None or dbpath is None:
			print("Both crystals directory and database path need to be specified.")
		else:
			crystals = os.listdir(crystal_dir)
			crystals = [os.path.normpath(crystal_dir+'/'+c) for c in crystals]
			#crystals = crystals[0:10000]

			mn = fastaai_consumer(dbpath)
			mn.open()
			mn.get_meta()
			mn.add_crystals(crystals)
			mn.close()

import argparse
def consumer_options():
	parser = argparse.ArgumentParser(description='Add genomes to a fastaai v2 database')
	
	parser.add_argument('-db', '--database', dest='database', default = None,
						help='')
	parser.add_argument('-c', '--crystals', dest = 'crystals_dir', default = None,
						help='Directory containing FastAAI v2 crystals as produced by fastaai2 preproc')


	args = parser.parse_known_args()
	return parser, vars(args[0])