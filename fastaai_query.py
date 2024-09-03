import sys
import os
import sqlite3
import numpy as np

from supports.metadata_loader import fastaai_metadata_loader

import multiprocessing
from multiprocessing.managers import SharedMemoryManager
from multiprocessing.shared_memory import SharedMemory

import ctypes

import datetime
import time

import json

import argparse

#Creates vectors organized as genome: acc: tetras
class fastaai_db_query:
	def __init__(self, dbpath):
		self.path = dbpath
		self.conn = None
		self.curs = None
		
		self.accessions = None
		
		self.metadata = None
		
		self.accession_index = None
		self.genome_index = None
		self.reverse_gix = None
		self.reverse_acc = None
		
		self.gak = None
		self.gas = None
		
		self.start_index = 0
		self.starts_by_genome = None
		self.ends_by_genome = None
		self.starts_within_genome = None
		self.ends_within_genome = None
		self.genome_record = None
		self.recovered_accessions = None
		
		self.genome_tax = None
		self.genlens = None
		
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

		#Query genomes are reclassified against the target.
		#self.gen_classes = self.metadata.genome_classes
		self.genlens = self.metadata.genlens
		self.genome_tax = self.metadata.gen_tax

		self.gak = self.metadata.gak #kmer counts
		self.gas = self.metadata.gas #hmm scores
		
		self.accession_index = self.metadata.scp_name_to_id
		self.genome_index = self.metadata.gix
		self.reverse_gix = self.metadata.reverse_gix

		self.reverse_acc = self.metadata.scp_id_to_name
		self.recovered_accessions = self.metadata.recovered_accessions
	
	def load_accessions_by_group(self, accessions = [], low_high_tups = None):
		self.starts_within_genome = {}
		self.ends_within_genome = {}

		self.genome_record = {}
		self.recovered_accessions = {}

		if low_high_tups is None:
			low_high_tups = [(0, len(self.genome_index),)]
		

		for lh in low_high_tups:
			low_bound = lh[0]
			high_bound = lh[1]
			
			lh_key = "query_"+str(low_bound)+"_to_"+str(high_bound)

			self.starts_within_genome[lh_key] = {}
			self.ends_within_genome[lh_key] = {}
			self.genome_record[lh_key] = []
			self.recovered_accessions[lh_key] = []

			self.load_accessions(accessions, low_bound, high_bound, lh_key)
			
	def load_accessions(self, acc_list = [], genome_range_low = 0, genome_range_high = None, group_key = None):
		self.start_index = 0
		
		recovered_data = {}
		
		for acc in acc_list:
			next_chunk = self.load_one_accession_genomes(acc, genome_range_low, genome_range_high, group_key)
			if next_chunk is not None:
				for g in next_chunk:
					if g not in recovered_data:
						recovered_data[g] = {}
					
					recovered_data[g][acc] = next_chunk[g]
					
		
		big_boi = []
		start_index = 0
		self.starts_within_genome[group_key] = {}
		self.ends_within_genome[group_key] = {}
		
		for genome in recovered_data:
			self.starts_within_genome[group_key][genome] = {}
			self.ends_within_genome[group_key][genome] = {}
			for acc in recovered_data[genome]:
				self.starts_within_genome[group_key][genome][acc] = start_index
				start_index += recovered_data[genome][acc].shape[0]
				big_boi.append(recovered_data[genome][acc])
				self.ends_within_genome[group_key][genome][acc] = start_index
				recovered_data[genome][acc] = None
				
			recovered_data[genome] = None
				
		recovered_data = None
		if len(big_boi) > 0:
			self.genome_record[group_key] = np.concatenate(big_boi)
			self.recovered_accessions[group_key] = set(self.recovered_accessions[group_key])
		else:
			self.genome_record[group_key] = None
			self.recovered_accessions[group_key] = None
		
	def load_one_accession_genomes(self, accession_name, low, high, group_key):
		sql = 'SELECT * FROM "{acc}_genomes" WHERE genome_id >= ? AND genome_id < ?'.format(acc=accession_name)
		results = self.curs.execute(sql, (low, high,)).fetchall()
		tetramer_record = None
		if len(results) > 0:
			tetramer_record = {}
			for r in results:
				genome_id = r[0]
				tetramer_list = np.frombuffer(r[1], dtype = np.int32)
				
				tetramer_record[genome_id] = tetramer_list
				
				self.recovered_accessions[group_key].append(accession_name)
			
		return tetramer_record
			
#Creates vectors of genomes organized as acc : tetra : genomes
class fastaai_db_target:
	def __init__(self, dbpath):
		self.path = dbpath
		self.conn = None
		self.curs = None
		
		self.accessions = None
		
		self.metadata = None
		
		self.accession_index = None
		self.genome_index = None
		self.reverse_gix = None
		self.reverse_acc = None
		
		self.gak = None
		self.gas = None
		
		self.start_index = 0
		self.starts_by_acc = None
		self.ends_by_acc = None
		self.starts_within_acc = None
		self.ends_within_acc = None
		self.genome_record = None
		self.recovered_accessions = None
		
		self.gen_classes = None
		self.genome_tax = None
		self.genlens = None
		
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
		#Genome-level data
		
		self.metadata = fastaai_metadata_loader(self.path)
		self.metadata.load_meta(skip_gak = False)

		self.gak = self.metadata.gak #kmer counts
		self.gas = self.metadata.gas #hmm scores
		
		self.gen_classes = self.metadata.genome_classes
		self.genlens = self.metadata.genlens
		self.genome_tax = self.metadata.gen_tax
		
		self.accession_index = self.metadata.scp_name_to_id
		self.genome_index = self.metadata.gix
		self.reverse_gix = self.metadata.reverse_gix

		self.reverse_acc = self.metadata.scp_id_to_name
		self.recovered_accessions = self.metadata.recovered_accessions
		
		self.genome_class_rules = None
		self.final_transform_rules = None
		
		self.aai_models = None
		self.aai_jacc_tf = None

		self.query_genome_classes = None

	def load_accessions_by_group(self, accessions, low_high_tups = None):
		self.starts_by_acc = {}
		self.ends_by_acc = {}
		self.starts_within_acc = {}
		self.ends_within_acc = {}
		self.genome_record = {}
		self.recovered_accessions = {}

		if low_high_tups is None:
			low_high_tups = [(0, len(self.genome_index),)]
		

		for lh in low_high_tups:
			low_bound = lh[0]
			high_bound = lh[1]
			
			lh_key = "target_"+str(low_bound)+"_to_"+str(high_bound)

			self.starts_by_acc[lh_key] = {}
			self.ends_by_acc[lh_key] = {}
			self.starts_within_acc[lh_key] = {}
			self.ends_within_acc[lh_key] = {}
			self.genome_record[lh_key] = []
			self.recovered_accessions[lh_key] = []

			self.load_accessions(accessions, low_bound, high_bound, lh_key)

	def load_accessions(self, acc_list = [], genome_range_low = 0, genome_range_high = None, group_key = None):
		self.start_index = 0
			
		for acc in acc_list:
			self.starts_by_acc[group_key][acc] = self.start_index
			#print("\tLoading:", acc)
			next_chunk, next_starts, next_ends = self.load_one_accession_genomes(acc, genome_range_low, genome_range_high, group_key)
			if next_chunk is not None:
				#switching from within to absolute values
				self.starts_within_acc[group_key][acc] = next_starts
				self.ends_within_acc[group_key][acc] = next_ends
				self.genome_record[group_key].append(next_chunk)
				self.start_index += next_chunk.shape[0]
			self.ends_by_acc[group_key][acc] = self.start_index
		
		self.genome_record[group_key] = np.concatenate(self.genome_record[group_key])
		
		self.recovered_accessions[group_key] = set(self.recovered_accessions[group_key])
					
	#Credit to Thomas Browne at 
	#https://stackoverflow.com/questions/1066758/find-length-of-sequences-of-identical-values-in-a-numpy-array-run-length-encodi	
	#returns tuple of run lengths, start positions, and values
	def rle(self, ia):
		""" run length encoding. Partial credit to R rle function. 
			Multi datatype arrays catered for including non Numpy
			returns: tuple (runlengths, startpositions, values) """
		#ia = np.asarray(inarray)                # force numpy - we're pre-checking this, so no need to recreate the array
		n = len(ia)
		if n == 0: 
			return (None, None, None)
		else:
			y = ia[1:] != ia[:-1]               # pairwise unequal (string safe)
			i = np.append(np.where(y), n - 1)   # must include last element posi
			z = np.diff(np.append(-1, i))       # run lengths
			p = np.cumsum(np.append(0, z))[:-1] # positions
			return(z, p, ia[i])
		
	def load_one_accession_genomes(self, accession_name, low, high, group_key):
		sql = 'SELECT * FROM "{acc}_genomes" WHERE genome_id >= ? AND genome_id < ?'.format(acc=accession_name)
		results = self.curs.execute(sql, (low, high,)).fetchall()
		per_genome_starts = {}
		per_genome_ends = {}
		per_tetra_starts = None
		per_tetra_ends = None
		local_index = 0
		if len(results) > 0:
			tetramer_record = []
			for r in results:
				genome_id = r[0]
				tetramer_list = np.frombuffer(r[1], dtype = np.int32)
				per_genome_starts[genome_id] = local_index
				local_index += tetramer_list.shape[0]
				per_genome_ends[genome_id] = local_index
				
				tetramer_record.append(tetramer_list)
			
			tetramer_record = np.concatenate(tetramer_record)
			if tetramer_record.size == 0:
				genome_record = None
			else:
				self.recovered_accessions[group_key].append(accession_name)
				#Push the sort to here so that the sort is implicitly within each accession and that logic stays
				#Create a vector of genome IDs
				genome_record = np.zeros(shape = tetramer_record.shape, dtype = np.int32)
				for genome_id in per_genome_starts:
					genome_record[per_genome_starts[genome_id]:per_genome_ends[genome_id]] = genome_id
					
				#Get the tetramer-based ordering of genome IDs, then sort both lists
				ordering = tetramer_record.argsort()
				genome_record = genome_record[ordering]
				tetramer_record = tetramer_record[ordering]
				
				#Now we need RLE of tetramers to give us local starts, ends for each tetra:
				rle_tuple = self.rle(tetramer_record)
				#rle_tuple = rle(tetramer_record)
				tetramer_record = None
				run_lengths = rle_tuple[0]
				start_positions = rle_tuple[1]
				tetramer_values = rle_tuple[2]
				
				#Well thats easy
				#end_positions = start_positions + run_lengths
				
				#We could probably make these numpy arrays instead.
				per_tetra_starts = {}
				per_tetra_ends = {}
				for s, l, t in zip(start_positions, run_lengths, tetramer_values):
					per_tetra_starts[t] = s + self.start_index #Absolute positions in the overall shared vector
					per_tetra_ends[t] = s+l + self.start_index #Absolute positions in the overall shared vector
				
		else:
			genome_record = None
		
		return genome_record, per_tetra_starts, per_tetra_ends

	#The target database determines all rule-based info.
	def load_genome_class_rules_and_final_tf(self, query_acc_idx):
		sql = "SELECT ruleset_json FROM rule_addition_fields"
		rule = self.curs.execute(sql).fetchone()[0]
		class_extraction = json.loads(rule)
		
		self.rule_rider = None
		self.genome_class_rules = {}
		self.final_transform_rules = {}
		
		if 'class_assignment' in class_extraction:
			if 'classes' in class_extraction['class_assignment']:
				for c in class_extraction['class_assignment']['classes']:
					class_id = class_extraction['class_assignment']['classes'][c]["class"]
					rationale = class_extraction['class_assignment']['classes'][c]["rationale"]
					
					og_size = len(class_extraction['class_assignment']['classes'][c]['members'])
					scp_membership = []
					for scp in class_extraction['class_assignment']['classes'][c]['members']:
						if scp in query_acc_idx:
							scp_membership.append(query_acc_idx[scp])
					
					scp_membership = set(scp_membership)
					
					self.genome_class_rules[c] = {'class_id':class_id, 'scp_members':scp_membership, 'rationale':rationale, 'original_size':og_size}
			
			#Has to have a rider or None
			self.rule_rider = class_extraction['class_assignment']['rider']

		if 'final_transform' in class_extraction:
			for genome_class_1 in class_extraction['final_transform']:
				self.final_transform_rules[int(genome_class_1)] = {}
				for genome_class_2 in class_extraction['final_transform'][genome_class_1]:
					
					this_dat = class_extraction['final_transform'][genome_class_1][genome_class_2]
					intercept = this_dat['intercept']
					jslope = this_dat['jacc_or_aai_slope']
					scp_ct_slope = this_dat['shared_scp_ct']
					scp_dev_slope = this_dat['scp_stddev']
					next_rule = (intercept, jslope,
								scp_ct_slope, scp_dev_slope,)
					
					self.final_transform_rules[int(genome_class_1)][int(genome_class_2)] = next_rule				


	#The transformation models always come from the target database
	#Load per-SCP jacc -> AAI models
	#Note: this assumes that the genome class numbers are identical across two databases. 
	#I should make a way of checking to ensure this is true, or simply re-classify incoming query genomes according to the target db
	def load_aai_models(self):
		self.aai_models = {}
		sql = "SELECT * FROM rulesets"
		recovered_scps = self.curs.execute(sql).fetchall()
		for row in recovered_scps:
			scp = row[0]
			scpid = self.accession_index[scp]
			c1, c2 = row[1], row[2] #genome classes
			lb, hb = row[3], row[4] #low bound, high bound
			weight = row[5]
			jacc_int, jacc_slope = row[6], row[7]
			minlen, lendiff = row[8], row[9]
			minscore, scorediff = row[10], row[11]
			
			if scpid not in self.aai_models:
				self.aai_models[scpid] = {}
			if c1 not in self.aai_models[scpid]:
				self.aai_models[scpid][c1] = {}
			if c2 not in self.aai_models[scpid][c1]:
				self.aai_models[scpid][c1][c2] = [[], [], []]

			
			next_rule = np.array([weight, jacc_int, jacc_slope, minlen, lendiff, minscore, scorediff], dtype = np.float_)	
			self.aai_models[scpid][c1][c2][0].append(lb)
			self.aai_models[scpid][c1][c2][1].append(hb)
			self.aai_models[scpid][c1][c2][2].append(next_rule)
			
		#here we need to prep these for jaccard numpy piecwise functions, which is 
		for scpid in self.aai_models:
			condition_set = []
			jaccard_slope_set = []	

			for c1 in self.aai_models[scpid]:
				for c2 in self.aai_models[scpid][c1]:
					self.aai_models[scpid][c1][c2][0] = np.array(self.aai_models[scpid][c1][c2][0], dtype = np.float_)
					self.aai_models[scpid][c1][c2][1] = np.array(self.aai_models[scpid][c1][c2][1], dtype = np.float_)
					self.aai_models[scpid][c1][c2][2] = np.vstack(self.aai_models[scpid][c1][c2][2], dtype = np.float_)
	
	def classify_query_genomes(self, query_gak):
		self.query_genome_classes = {}
		for genome in query_gak:
			query_scps = set(query_gak[genome].keys())
			winning_membership = 0
			winning_class = 0
			for genome_class in self.genome_class_rules:
				membership_pct = len(query_scps.intersection(self.genome_class_rules[genome_class]['scp_members'])) / self.genome_class_rules[genome_class]['original_size']
				if membership_pct > winning_membership:
					winning_membership = membership_pct
					winning_class = self.genome_class_rules[genome_class]['class_id']
			
			self.query_genome_classes[genome] = winning_class
					

def split_indicies(max_val, num_grps):
	starts, ends = [], []
	paired_se = []
	splitsize = 1.0/num_grps*max_val
	for i in range(num_grps):
		#starts.append(round(i*splitsize))
		#ends.append(round((i+1)*splitsize))
		next_group = (int(round(i*splitsize)), int(round((i+1)*splitsize)), )
		paired_se.append(next_group)

	return paired_se

#Function to collect metadata for calculating Jaccard indices
#Needs to collect genome classes, lengths, HMM scores
def prep_target_genome_metadata(db_obj, accessions, starts_and_ends):
	sz = len(db_obj.genome_index)
	arrs = {}
	gas_arrs = {}
	
	for scp in accessions:
		scpid = db_obj.accession_index[scp]
		arrs[scpid] = np.zeros(sz, dtype = np.int32)
		gas_arrs[scpid] = np.zeros(sz, dtype = np.float_)
	for genomeid in db_obj.gak:
		for scpid in db_obj.gak[genomeid]:
			arrs[scpid][genomeid] = db_obj.gak[genomeid][scpid]
			gas_arrs[scpid][genomeid] = db_obj.gas[genomeid][scpid]
	
	recovered_scps = []
	numpyized_gak = []
	numpyized_gas = []
	for scpid in sorted(arrs):
		numpyized_gak.append(arrs[scpid])
		numpyized_gas.append(gas_arrs[scpid])
		
		
	recovered_scps = np.array(recovered_scps, dtype = np.int32)
	numpyized_gak = np.vstack(numpyized_gak)
	numpyized_gak = numpyized_gak.astype(np.int32)
	numpyized_gas = np.vstack(numpyized_gas)
	numpyized_gas = numpyized_gas.astype(np.float_)
	
	presabs = numpyized_gak > 0
	
	return numpyized_gak, presabs, numpyized_gas


def check_shared(names):
	qname = names[0]
	tname = names[1]
	target_start = names[2]
	target_end = names[3]
	
	print("Starting block", qname, "vs", tname, file = sys.stderr)

	sm = SharedMemory(qname)
	arr = np.frombuffer(sm.buf, dtype = np.int32)
	query_data = np.copy(arr)
	arr = None
	sm.close()

	sm = SharedMemory(tname)
	arr = np.frombuffer(sm.buf, dtype = np.int32)
	target_data = np.copy(arr)
	arr = None
	sm.close()
	
	sm = SharedMemory('target_gak')
	arr = np.frombuffer(sm.buf, dtype = np.int32)
	target_gak = np.copy(arr)
	arr = None
	sm.close()
	
	sm = SharedMemory('target_gas')
	arr = np.frombuffer(sm.buf, dtype = np.float_)
	target_gas = np.copy(arr)
	arr = None
	sm.close()
	
	sm = SharedMemory('target_presabs')
	arr = np.frombuffer(sm.buf, dtype = bool)
	target_presabs = np.copy(arr)
	arr = None
	sm.close()
	
	sm = SharedMemory('target_genclasses')
	arr = np.frombuffer(sm.buf, dtype = np.int32)
	target_genclass = np.copy(arr)
	arr = None
	sm.close()
	
	sm = SharedMemory('target_genlens')
	arr = np.frombuffer(sm.buf, dtype = np.int32)
	target_genlen = np.copy(arr)
	arr = None
	sm.close()
	
	#frombuffer is always 1d arr, reshape to match og.
	
	target_gak = target_gak.reshape(shape_me)
	target_gas = target_gas.reshape(shape_me)
	target_presabs = target_presabs.reshape(shape_me)
	
	#Subset to relevant target genomes
	target_gak = target_gak[:, target_start:target_end]
	target_gas = target_gas[:, target_start:target_end]
	target_presabs = target_presabs[:, target_start:target_end]
	
	target_genlen = target_genlen[target_start:target_end]
	target_genclass = target_genclass[target_start:target_end]
		
	unqiue_tgenclass = np.unique(target_genclass)
	uidx = {}
	for idx in unqiue_tgenclass:
		uidx[idx] = np.where(target_genclass == idx)[0]
		
	#query is vectors organized as genome: acc: tetras
	
	#target is vectors organized as acc : tetra : genomes

	starts_within_genome, ends_within_genome = qstarts[qname], qends[qname]
	starts_within_acc, ends_within_acc = tstarts[tname], tends[tname]
	
	outwriter = open(os.path.normpath("fastaai_output/"+qname+"_vs_"+tname+".txt"), "w")
	
	for genome in sorted(starts_within_genome):
		tgak_slice = []
		tpres_slice = []
		tgas_slice = []
		this_genome = []
		this_genome_unions = []
		
		query_genome_length = query_genlens[genome]
		query_tax = query_genome_taxonomy[genome]
		query_gclass = query_genclasses[genome]
		
		q_acc_scores = []
		
		observed_scps = []
		
		for acc in starts_within_genome[genome]:
			loaded_data = []
			bincts = None
			
			q_acc_id = query_accession_index[acc]
			t_acc_id = idx_translator[q_acc_id]
						
			if acc in starts_within_acc:
				startidx_query = starts_within_genome[genome][acc]
				endidx_query = ends_within_genome[genome][acc]
				
				for tetramer in query_data[startidx_query:endidx_query]:
					if tetramer in starts_within_acc[acc]:
						startidx_target = starts_within_acc[acc][tetramer]
						endidx_target = ends_within_acc[acc][tetramer]
						data = target_data[startidx_target:endidx_target]
						loaded_data.append(data)
			
			if len(loaded_data) > 0:
				q_acc_scores.append(query_acc_hmmscores[genome][q_acc_id])
				observed_scps.append(t_acc_id)
				
				loaded_data = np.concatenate(loaded_data)
				bincts = np.bincount(loaded_data, minlength = target_end + 1)
				bincts = bincts[target_start:target_end]
				this_genome.append(bincts)
				#Take target kmer counts and add query kmer count per SCP
				tgak_slice.append(target_gak[t_acc_id] + endidx_query - startidx_query)
				#tgak_slice.append(target_gak[t_acc_id] + endidx_query - startidx_query)
				tpres_slice.append(target_presabs[t_acc_id])
				tgas_slice.append(target_gas[t_acc_id])
				
		
		if len(this_genome) > 0:
			this_genome = np.vstack(this_genome)
			tgak_slice = np.vstack(tgak_slice)
			tpres_slice = np.vstack(tpres_slice)
			tgas_slice = np.vstack(tgas_slice)
			q_acc_scores = np.array(q_acc_scores, dtype = np.float_)
			
			shared_scp_ct = np.sum(tpres_slice, axis = 0)
			
			unions = tgak_slice - this_genome
			jaccards = this_genome / unions
			

			#here is where we have to apply the per-jaccard models.
			min_hmm_score = np.minimum(q_acc_scores[:, None], tgas_slice)
			max_hmm_score = np.maximum(q_acc_scores[:, None], tgas_slice)
			
			score_diff = max_hmm_score - min_hmm_score
			
			
			min_genome_length = np.minimum(query_genome_length, target_genlen)
			max_genome_length = np.maximum(query_genome_length, target_genlen)
			genome_length_diff = max_genome_length - min_genome_length
			

			#jacc_sums = np.sum(jaccards, axis = 0)
			#avg_jacc = jacc_sums / shared_scp_ct
			converted_aai = []
			converted_weights = []
			
			for i, scp in zip(range(0, len(observed_scps)), observed_scps):
				aai_model = aai_models[scp]
				calculated_jaccs = jaccards[i]
				transformed_aai = np.zeros(calculated_jaccs.shape, dtype = np.float_)
				weights_vector = np.zeros(calculated_jaccs.shape, dtype = np.float_)
				
				minscore_row = min_hmm_score[i]
				scorediff_row = score_diff[i]
				
				for idx in uidx:
					relevant_jaccs = calculated_jaccs[uidx[idx]]
					if query_gclass in aai_model:
						if idx in aai_model[query_gclass]:
							piecewise_linmods = aai_model[query_gclass][idx][2]
							#These are constant across the piecewise models, so row 0 is always OK
							minlength_mult = piecewise_linmods[0, 3]
							lengthdiff_mult = piecewise_linmods[0, 4]
							minscore_mult = piecewise_linmods[0, 5]
							scorediff_mult = piecewise_linmods[0, 6]
							
							minlength_addition = minlength_mult * min_genome_length[uidx[idx]]
							lengthdiff_addition = lengthdiff_mult * genome_length_diff[uidx[idx]]
							minscore_addition = minscore_mult * minscore_row[uidx[idx]]
							scorediff_addition = scorediff_mult * scorediff_row[uidx[idx]]
							
							jacc_slope_value = piecewise_linmods[:, 2]
							intercept_values = piecewise_linmods[:, 1]
							weights = piecewise_linmods[:, 0]
							#lowbound = aai_model[query_gclass][idx][0]
							highbound = aai_model[query_gclass][idx][1]
							inserts = np.searchsorted(highbound, relevant_jaccs, side = 'left')
							
							#add piecewise jaccard slope * observed jaccard
							transformed_aai[uidx[idx]] = np.multiply(relevant_jaccs, jacc_slope_value[inserts])
							#add intercept
							transformed_aai[uidx[idx]] += intercept_values[inserts]
							
							#Add piecewise local weight
							weights_vector[uidx[idx]] += weights[inserts]
							#add other slope info
							transformed_aai[uidx[idx]] += (minlength_addition + lengthdiff_addition + minscore_addition + scorediff_addition)
						
				converted_aai.append(transformed_aai)
				converted_weights.append(weights_vector)
				
			converted_aai = np.vstack(converted_aai)
			converted_weights = np.vstack(converted_weights)
			
			aai_predictions = np.multiply(converted_aai, converted_weights)
			aai_predictions = np.sum(aai_predictions, axis = 0)
			
			weight_sums = np.sum(converted_weights, axis = 0)
			aai_predictions = aai_predictions / weight_sums

			#standard deviations taken over only the values where an intersection was found
			aai_deviation = np.std(converted_aai, axis = 0, where = tpres_slice)
			
			np.nan_to_num(aai_predictions, copy = False)
			np.nan_to_num(aai_deviation, copy = False)

			#Final model adjustments
			if query_gclass in final_tf:
				for idx in uidx:
					if idx in final_tf[query_gclass]:
						#Final rules like so
						#next_rule = (intercept, jslope,
						#scp_ct_slope, scp_dev_slope) 
						
						model = final_tf[query_gclass][idx]
						relevant_aai = aai_predictions[uidx[idx]]
						relevant_stdev = aai_deviation[uidx[idx]]
						relevant_scp_cts = shared_scp_ct[uidx[idx]]
						
						#print(query_gclass, idx, relevant_aai, final_tf[query_gclass][idx])
						
						
						relevant_aai = (final_tf[query_gclass][idx][0] +  #intercept
						(relevant_aai * final_tf[query_gclass][idx][1]) +  #aai_est slope
						
						(relevant_scp_cts * final_tf[query_gclass][idx][2]) + #shared scp_slope
						(relevant_stdev * final_tf[query_gclass][idx][3])) #stddev slope
						
						#print(relevant_aai)
						
						aai_predictions[uidx[idx]] = relevant_aai
						
			aai_predictions[aai_predictions > 100.0] = 100.0 #correct to make sense.
			
			query_genome_name = query_genidx[genome]
			for i, tgt in zip(range(0, len(aai_predictions)), range(target_start, target_end)):
				target_genome_name = target_genidx[tgt]
				target_tax = target_genome_taxonomy[tgt]
				print(query_genome_name, query_tax, target_genome_name, target_tax, aai_predictions[i], aai_deviation[i], shared_scp_ct[i], sep = "\t", file = outwriter)
			
				
		else:
			query_genome_name = query_genidx[genome]
			for i, tgt in zip(range(0, len(aai_predictions)), range(target_start, target_end)):
				target_genome_name = target_genidx(tgt)
				print(query_genome_name, query_tax, target_genome_name, target_tax, 0, 0, 0, sep = "\t", file = outwriter)

	outwriter.close()
	
	return (qname, tname,)


def run_fastaai_query():
	print("Program start at", datetime.datetime.now(), file = sys.stderr)
		
	parser, args = query_exec_options()
	
	if len(sys.argv) < 3:
		print(parser.print_help())
	else:
		
		qdb = args['query_database']
		tdb = args['target_database']
		num_thds = args['threads']
		query_chunk_size = args['query_chunk_size']
		tgt_chunk_size = args['target_chunk_size']

		qmn = fastaai_db_query(qdb)
		tmn = fastaai_db_target(tdb)
			
		print("Loading query database metadata", file = sys.stderr)
		qmn.open()
		qmn.get_meta()
		qmn.close()

		print("Loading target database metadata", file = sys.stderr)
		tmn.open()
		tmn.get_meta()
		tmn.load_genome_class_rules_and_final_tf(qmn.accession_index)
		tmn.load_aai_models()
		tmn.classify_query_genomes(qmn.gak)
		tmn.close()

		shared_scps = qmn.recovered_accessions.intersection(tmn.recovered_accessions)
		shared_scps = list(shared_scps)
		shared_scps.sort()

		num_queries = len(qmn.genome_index)
		num_tgts = len(tmn.genome_index)
		
		#Default to load a reasonable chunking
		if num_thds > 1:
			if tgt_chunk_size == 0:
				tgt_chunk_size = int(num_queries / int(num_thds / 2))
			if query_chunk_size == 0:
				query_chunk_size = int(num_tgts / (num_thds - int(num_thds / 2)))
		else:
			if tgt_chunk_size == 0:
				tgt_chunk_size = num_tgts
			if query_chunk_size == 0:
				query_chunk_size = num_queries



		#Ensure chunks are possible
		if query_chunk_size > num_queries:
			print("Fewer query genomes in the database ({sz1}) than the requested query chunk size ({sz2}).".format(sz1 = str(num_queries), sz2 = str(query_chunk_size)), file = sys.stderr)
			print("Reducing query chunk size to number of database genomes.", file = sys.stderr)
			query_chunk_size = num_queries
			
		if tgt_chunk_size > num_tgts:
			print("Fewer target genomes in the database ({sz1}) than the requested target chunk size ({sz2}).".format(sz1 = str(num_tgts), sz2 = str(tgt_chunk_size)), file = sys.stderr)
			print("Reducing target chunk size to number of database genomes.", file = sys.stderr)
			tgt_chunk_size = num_tgts
		
		#print(query_chunk_size)
		#print(tgt_chunk_size)
		#quit()

		start_end_pairs_queries = split_indicies(num_queries, int(num_queries/query_chunk_size))
		start_end_pairs_targets = split_indicies(num_tgts, int(num_tgts/tgt_chunk_size))

		print("Loading query database", file = sys.stderr)
		qmn.open()
		qmn.load_accessions_by_group(shared_scps, start_end_pairs_queries)
		qmn.close()

		print("Loading target database", file = sys.stderr)
		tmn.open()
		tmn.load_accessions_by_group(shared_scps, start_end_pairs_targets)
		tmn.close()

		#global qgak
		#global tgak
		#global qpresabs
		#global t_accidx

		global query_genidx
		global target_genidx

		query_genidx = qmn.metadata.reverse_gix
		target_genidx = tmn.metadata.reverse_gix

		global query_accession_index
		global idx_translator
		global shape_me

		#global tpresabs

		global qstarts
		global qends
		global tstarts
		global tends


		global query_genclasses
		global query_genlens
		global query_acc_hmmscores

		global query_genome_taxonomy

		#global target_genclasses
		#global target_genlens
		#global target_acc_hmmscores
		global target_genome_taxonomy

		#Last bits
		global aai_models
		global final_tf

		#We use the target database to classify incoming query genomes to ensure they match target database expectations 
		query_genclasses = tmn.query_genome_classes
		#And now back to the query database
		query_genlens = qmn.genlens
		query_acc_hmmscores = qmn.gas
		query_genome_taxonomy = qmn.genome_tax

		query_accession_index = qmn.accession_index

		#We're gonna want to convert the target ones to be vectors
		target_genclasses = tmn.gen_classes
		target_genlens = tmn.genlens

		#This one can stay as a dict
		target_genome_taxonomy = tmn.genome_tax

		#both tax dicts should be made into gtdb fmt strings...
		def gtdb_fmt_tax(tax_dict):
			formatted_tax = {}
			for genome in tax_dict:
				next_tax = []
				for rank, val in zip(["d__", "p__", "c__", "o__", "f__", "g__", "s__"], tax_dict[genome]):
					if val is not None:
						next_tax.append(rank + val)
					else:
						next_tax.append(rank+"")
				next_tax = ';'.join(next_tax)
				formatted_tax[genome] = next_tax
			
			return formatted_tax
			
		query_genome_taxonomy = gtdb_fmt_tax(query_genome_taxonomy)
		target_genome_taxonomy = gtdb_fmt_tax(target_genome_taxonomy)


		tgclass = []
		for t in sorted(target_genclasses.keys()):
			tgclass.append(target_genclasses[t])
			
		tgclass = np.array(tgclass, dtype = np.int32)
		target_genclasses = tgclass

		tglens = []
		for t in sorted(target_genlens.keys()):
			tglens.append(target_genlens[t])
			
		tglens = np.array(tglens, dtype = np.int32)
		target_genlens = tglens

		aai_models = tmn.aai_models
		final_tf = tmn.final_transform_rules

		#We search query -> target, so we implicitly have query data every time
		#qgak, qpresabs = prep_genome_metadata(qmn, shared_scps, start_end_pairs_queries)

		#target data we're gonna put into shared memory.

		idx_translator = {}
		for scp in shared_scps:
			qidx = qmn.accession_index[scp]
			tidx = tmn.accession_index[scp]
			idx_translator[qidx] = tidx

		tgak, tpresabs, tgas = prep_target_genome_metadata(tmn, shared_scps, start_end_pairs_targets)
		shape_me = tgak.shape

		qstarts = qmn.starts_within_genome
		qends = qmn.ends_within_genome

		tstarts = tmn.starts_within_acc
		tends = tmn.ends_within_acc


		args = []
		for qname in qmn.genome_record:
			for tname, se in zip(tmn.genome_record, start_end_pairs_targets):
				next_pair = (qname, tname, se[0], se[1],)
				args.append(next_pair)

		if not os.path.exists("fastaai_output"):
			os.mkdir("fastaai_output")

		print("Starting calculation at", datetime.datetime.now(), file = sys.stderr)
		with SharedMemoryManager() as memory_manager:
			shared_memories = []
			memory_names = []
			
			shared_memory = SharedMemory(name = 'target_gak', create=True, size = tgak.nbytes)
			data = np.ndarray(tgak.shape, dtype = np.int32, buffer = shared_memory.buf)
			data[:] = tgak
			shared_memories.append(shared_memory)
			memory_names.append('target_gak')
			
			shared_memory = SharedMemory(name = 'target_presabs', create=True, size = tpresabs.nbytes)
			data = np.ndarray(tpresabs.shape, dtype = bool, buffer = shared_memory.buf)
			data[:] = tpresabs
			shared_memories.append(shared_memory)
			memory_names.append('target_presabs')
			
			shared_memory = SharedMemory(name = 'target_gas', create=True, size = tgas.nbytes)
			data = np.ndarray(tgas.shape, dtype = np.float_, buffer = shared_memory.buf)
			data[:] = tgas
			shared_memories.append(shared_memory)
			memory_names.append('target_gas')
			
			shared_memory = SharedMemory(name = 'target_genlens', create=True, size = target_genlens.nbytes)
			data = np.ndarray(target_genlens.shape, dtype = np.int32, buffer = shared_memory.buf)
			data[:] = target_genlens
			shared_memories.append(shared_memory)
			memory_names.append('target_genlens')
			
			shared_memory = SharedMemory(name = 'target_genclasses', create=True, size = target_genclasses.nbytes)
			data = np.ndarray(target_genclasses.shape, dtype = np.int32, buffer = shared_memory.buf)
			data[:] = target_genclasses
			shared_memories.append(shared_memory)
			memory_names.append('target_genclasses')

			#Load query data to shared mem
			for qname in qmn.genome_record:
				shared_memory = SharedMemory(name = qname, create=True, size = qmn.genome_record[qname].nbytes)
				data = np.ndarray(qmn.genome_record[qname].shape, dtype = np.int32, buffer = shared_memory.buf)
				data[:] = qmn.genome_record[qname]
				shared_memories.append(shared_memory)
				memory_names.append(qname)
				qmn.genome_record[qname] = None
			
			#load target data to shared mem
			for tname in tmn.genome_record:
				shared_memory = SharedMemory(name = tname, create=True, size = tmn.genome_record[tname].nbytes)
				data = np.ndarray(tmn.genome_record[tname].shape, dtype = np.int32, buffer = shared_memory.buf)
				data[:] = tmn.genome_record[tname]
				shared_memories.append(shared_memory)
				memory_names.append(tname)
				tmn.genome_record[tname] = None

			pool = multiprocessing.Pool(num_thds, maxtasksperchild = 1)
			for result in pool.imap_unordered(check_shared, args):
				print(result, "complete", file = sys.stderr)

			pool.close()
			pool.join()

			#clean up or mem leak
			for s in shared_memories:
				s.close()
				s.unlink()

			shared_memories = None
			memory_names = None



def query_exec_options():
	#needs
	#qdb = args['query_database']
	#tdb = args['target_database']
	#num_thds = args['threads']
	#query_chunk_size = args['query_chunk_size']
	#tgt_chunk_size = args['target_chunk_size']
	
	parser = argparse.ArgumentParser(description='Options for executing a FastAAI v2 database search. --query_database and --target_database are required')
	
	parser.add_argument('-qdb', '--query_database', dest = 'query_database',
						help='')
	parser.add_argument('-tdb', '--target_database', dest='target_database',
						help='')
	parser.add_argument('-t', '--threads', dest='threads', default=1, type = int,
						help='Number of processes to use.')
	parser.add_argument('--query_chunk_size', dest='query_chunk_size', default = 0, type = int,
						help='The number of query genomes to load at a time. Larger chunk size is faster, but consumes more memory.')
						
	parser.add_argument('--target_chunk_size', dest='target_chunk_size', default = 0, type = int,
						help='The number of target genomes to load at a time. Larger chunk size is faster, but consumes more memory. Larger target chunk size is more efficient than larger query size.')

	args = parser.parse_known_args()
	return parser, vars(args[0])

