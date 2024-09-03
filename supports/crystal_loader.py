import sys
import os
import numpy as np
import json

try:
	from .agnostic_reader import agnostic_reader
	from .get_file_basename import get_basename
except:
	from agnostic_reader import agnostic_reader
	from get_file_basename import get_basename

#Class for loading crystal files into JSON objects to be consumed by a FastAAI database
class ravenous_crystal_lizard:
	def __init__(self):
		self.db_genome_index = {}
		self.db_hmm_index = {}
		self.next_gen_id = 0
	
		self.gen_desc = None
		self.json_rep = None
		self.per_scp_tetras = None
		
		self.genome_name = None
		self.genome_length = None
		self.prot_basect = None
		self.num_prots = None
		
		self.taxonomy = None
		self.tax_domain = None
		self.tax_phylum = None
		self.tax_class = None
		self.tax_order = None
		self.tax_family = None
		self.tax_genus = None
		self.tax_species = None
		
	def set_indices(self, genome_index, hmm_index, next_gen_id):
		self.db_genome_index = genome_index
		self.db_hmm_index = hmm_index
		self.next_gen_id = next_gen_id
		
	def digest_taxonomy(self):
		tax_levels = self.taxonomy.split(";") #Split GTDB string
		tax_levels = [t.split("__")[1] for t in tax_levels] #Split each tax level to remove the tax rank prefix
		for i in range(0, len(tax_levels)):
			if len(tax_levels[i]) == 0:
				tax_levels[i] = None
		
		self.tax_domain  = tax_levels[0]
		self.tax_phylum  = tax_levels[1]
		self.tax_class   = tax_levels[2]
		self.tax_order   = tax_levels[3]
		self.tax_family  = tax_levels[4]
		self.tax_genus   = tax_levels[5]
		self.tax_species = tax_levels[6]
		
		
	def load_crystal(self, file):
		ar = agnostic_reader(file)
		dat = ar.read()
		ar = None
		as_json = json.loads(dat)
		self.genome_name = as_json["genome"]
		
		if 'genome_length' in as_json:
			self.genome_length = as_json["genome_length"]
		else:
			self.genome_length =  -1
		
		if 'coding_bases' in as_json:
			self.prot_basect = as_json['coding_bases']
		else:
			self.prot_basect = -1
			
		if 'protein_count' in as_json:
			self.num_prots = as_json['protein_count']
		else:
			self.num_prots = len(as_json["SCPs"])
			
		if 'taxonomy' in as_json:
			self.taxonomy = as_json['taxonomy']
		else:
			self.taxonomy = 'd__;p__;c__;o__;f__;g__;s__'
			
		self.digest_taxonomy()
	
		
		if self.genome_name not in self.db_genome_index:
			self.db_genome_index[self.genome_name] = self.next_gen_id
			self.next_gen_id += 1
			
		genome_id = self.db_genome_index[self.genome_name]
		
		self.gen_desc = []
		
		'''
		#CREATE TABLE genome_metadata 
			(genome_name TEXT, 
			genome_id INTEGER PRIMARY KEY, 
			genome_length INTEGER, 
			proteome_ct INTEGER, 
			scp_count INTEGER, 
			genome_class INTEGER,
			genome_class_label TEXT,
			dom TEXT, phyl TEXT, cl TEXT, ord TEXT, fam TEXT, genu TEXT, spec TEXT)
		'''

		
		self.json_rep = []
		
		self.per_scp_tetras = {}
		
		for protein in as_json["SCPs"]:
			'''
			per-protein JSON data as:
			"hmm_accession":self.protein_to_accs[p],
			"hmm_score":self.accs_to_besthits[self.protein_to_accs[p]],
			"protein_len":self.seqlens[p],
			"valid_tetramer_ct":self.kmer_cts[p],
			"base_20_tetras":self.proteome_seq[p]
			'''
			tetramers = np.array(as_json["SCPs"][protein]['base_20_tetras'], dtype = np.int32).tobytes()
			score = np.float_(as_json["SCPs"][protein]['hmm_score'])
			acc = as_json["SCPs"][protein]['hmm_accession']
			plen = as_json["SCPs"][protein]['protein_len']
			tct = as_json["SCPs"][protein]['valid_tetramer_ct']
			
			#if acc not in self.db_hmm_index:
				#print("HMM Accession", acc, "not found within the database. This protein will be skipped.")
				
			if acc in self.db_hmm_index:
				acc_id = self.db_hmm_index[acc]
			
				'''
				#this is the primary storage for the protein genome:accession:kmers (GAK) table
				
				"genome_id INTEGER",
				"scp_id INTEGER",
				"protein_length INTEGER",
				"hmm_score REAL",
				"tetra_count INTEGER"
				'''
				next_gak = [genome_id, acc_id, plen, score, tct]
				
				next_scp = [genome_id, tetramers]
				
				self.per_scp_tetras[acc] = next_scp
				self.json_rep.append(next_gak)
		
		if self.genome_length > 0:
			genlen_to_add = self.genome_length
		else:
			genlen_to_add = self.prot_basect
			
		#This is how the database wants the genome's info
		genome_metadata_addition = [self.genome_name,
									genome_id,
									genlen_to_add, #a minimum possible genome length based on number of bases found in all proteins; should be close for prokaryotes
									self.num_prots,
									len(self.json_rep), #VALID SCPs
									0, 'default_genome_class', #genome classes
									self.tax_domain, self.tax_phylum, self.tax_class, self.tax_order, self.tax_family, self.tax_genus, self.tax_species] #taxonomy placeholders
		
		self.gen_desc = genome_metadata_addition
		
		#Corresponds to tables
		#genome_metadata, genome_accession_counts, {scp}_genome
		return self.gen_desc, self.json_rep, self.per_scp_tetras
		
