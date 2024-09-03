import sys
import os

import pyhmmer

import gzip
import io

from .agnostic_reader import agnostic_reader
from .fasta_importer import fasta_file

from .get_file_basename import get_basename

import numpy as np

import multiprocessing

import sqlite3

class new_pyhmmer_manager:
	def __init__(self, compress = False):
		self.proteins = None
		self.hmm_models = None
		
		self.nt_easel = pyhmmer.easel.Alphabet.dna()
		self.amino_easel = pyhmmer.easel.Alphabet.amino()
		
		self.proteins_to_search = None
		
		self.hmm_result_proteins = []
		self.hmm_result_accessions = []
		self.hmm_result_scores = []
		
		self.printable_lines = []
		
		self.best_hits = None
		
		self.do_compress = compress
		
		self.pyhmmer_text = None
		
		self.acc_to_score = None
		
	def load_hmm_from_file(self, hmm_file):
		ar = agnostic_reader(hmm_file)
		hmm_text = ar.read()
		hmm_io = io.BytesIO(hmm_text.encode(encoding = "ascii"))
		self.convert_hmm_to_digital(hmm_io)
			
	def load_hmm_from_database(self, database):
		if os.path.exists(database):
			conn = sqlite3.connect(database)
			curs = conn.cursor()
			
			hmm_text = []
			sql = "SELECT scp_model FROM scp_data"
			for hmm in curs.execute(sql).fetchall():
				next_model = hmm[0]
				next_model = next_model.replace('"', '') #SQL version encloses the text with quotations so it's always a string
				next_model = next_model.replace('\\n', '\n') #SQL treats newlines as literal characters via double escape.
				hmm_text.append(next_model)
			
			curs.close()
			conn.close()
			
			hmm_text = '\n'.join(hmm_text)
			
			hmm_io = io.BytesIO(hmm_text.encode(encoding = "ascii"))
						
			self.convert_hmm_to_digital(hmm_io)
			
		else:
			print("Database", database, "not found. Cannot search prokaryotic HMMs")
		
	def convert_hmm_to_digital(self, hmm_io):
		with pyhmmer.plan7.HMMFile(hmm_io) as fh:
			hmms = list(fh)
			
		self.hmm_models = hmms
		#for h in self.hmm_models:
		#	print(h.accession)
		
	def load_protein_seqs_from_file(self, protein_file):
		ff = fasta_file(protein_file)
		ff.import_fasta()
		
		seqdict = ff.sequences
		
		self.convert_protein_seqs_in_mem(seqdict)
		
	def convert_protein_seqs_in_mem(self, contents):
		#Clean up.
		self.proteins_to_search = []
		
		for protein in contents:
			#Skip a protein if it's longer than 100k AA.
			if len(contents[protein]) >= 100000:
				continue
			as_bytes = protein.encode()
			#Pyhmmer digitization of sequences for searching.
			easel_seq = pyhmmer.easel.TextSequence(name = as_bytes, sequence = contents[protein])
			easel_seq = easel_seq.digitize(self.amino_easel)
			self.proteins_to_search.append(easel_seq)
			
		easel_seq = None
		
	def search_protein(self):
		top_hits = list(pyhmmer.hmmsearch(self.hmm_models, self.proteins_to_search, cpus=1, bit_cutoffs="trusted"))
		
		self.pyhmmer_text = io.BytesIO(b"")
		is_first = True
		for model in top_hits:
			if is_first:
				model.write(self.pyhmmer_text, header = True)
				is_first = False
			else:
				model.write(self.pyhmmer_text, header = False)
		

		#self.printable_lines = []
		
		self.hmm_result_proteins = []
		self.hmm_result_accessions = []
		self.hmm_result_scores = []
		
		for model in top_hits:
			for hit in model:
				target_name = hit.name.decode()
				query_acc = hit.best_domain.alignment.hmm_accession.decode()
				best_dom_score = hit.best_domain.alignment.domain.score

				self.hmm_result_proteins.append(target_name)
				self.hmm_result_accessions.append(query_acc)
				self.hmm_result_scores.append(best_dom_score)
						
	def filter_to_best_hits(self, best_hit_per_scp = True, reciprocal_best_hit = True):
		hmm_file = np.transpose(np.array([self.hmm_result_proteins, self.hmm_result_accessions, self.hmm_result_scores]))
				
		#hmm_file = np.loadtxt(hmm_file_name, comments = '#', usecols = (0, 3, 8), dtype=(str))
		#Sort the hmm file based on the score column in descending order.
		hmm_file = hmm_file[hmm_file[:,2].astype(float).argsort()[::-1]]
		
		if best_hit_per_scp:
			#Filter the file again for the unique ACCESSION names, since we're only allowed one gene per accession, I guess?
			#Don't sort the indices, we don't care about the scores anymore.
			hmm_file = hmm_file[np.sort(np.unique(hmm_file[:,1], return_index = True)[1])]
			
		if reciprocal_best_hit:
			#Identify the first row where each gene name appears, after sorting by score; 
			#in effect, return the highest scoring assignment per gene name
			#Sort the indices of the result to match the score-sorted table instead of alphabetical order of gene names
			hmm_file = hmm_file[np.sort(np.unique(hmm_file[:,0], return_index = True)[1])]
		
		#sql_friendly_names = [i.replace(".", "_") for i in hmm_file[:,1]]
		#self.best_hits = dict(zip(hmm_file[:,0], sql_friendly_names))
		self.best_hits = dict(zip(hmm_file[:,0], hmm_file[:,1]))

		self.acc_to_score = dict(zip(hmm_file[:,1], hmm_file[:,2].astype(float)))
		hmm_file = None
		
		return self.best_hits, self.acc_to_score

	def to_hmm_file_pyhmmer_method(self, output):
		#PyHMMER data is a bit hard to parse. For each result:
		content = self.pyhmmer_text.getvalue()

		if self.do_compress:						
			fh = gzip.open(output+".gz", "wb")
			fh.write(content)
			fh.close()
			content = None
			
		else:
			content = content.decode(encoding = "ascii")
			fh = open(output, "w")
			fh.write(content)
			fh.close()
			
		content = None

	


	
