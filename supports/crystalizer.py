import sys
import os
import numpy as np
import json
import multiprocessing

from .fasta_importer import fasta_file
from .protein_to_tetras import tetra_mapper
from .hmm_results_reader import hmm_results

from .get_file_basename import get_basename

import gzip

class crystalizer:
	def __init__(self, taxonomy = 'd__NA;p__NA;c__NA;o__NA;f__NA;g__NA;s__NA', compress = False):
		self.pf = None
		self.hf = None
		self.gf = None
		
		self.basename = None
		
		self.genome_length = None
		
		self.protein_to_accs = None
		self.accs_to_besthits = None
		
		self.proteome_seq = None
		
		self.tt = tetra_mapper()
		self.seqlens = {}
		self.kmer_cts = {}
		
		self.json = None
		
		self.min_genome_length = 0
		self.protein_count = 0
		
		self.taxonomy = taxonomy
		
		self.compress = compress
		
	def load_genome(self):
		#Default case
		self.genome_length = 0
		if self.gf is not None:
			mn = fasta_file(self.gf)
			mn.import_fasta()
			for s in mn.sequences:
				self.genome_length += len(mn.sequences[s])
		mn = None
		
	def load_hmm(self):
		mn = hmm_results(self.hf)
		self.protein_to_accs, self.accs_to_besthits = mn.run()
		
	def load_prots(self):
		mn = fasta_file(self.pf)
		mn.import_fasta()
		self.proteome_seq = mn.sequences
		
		mn = None
			
		
	def filter_prots(self):
		prots = list(self.proteome_seq.keys())
		self.protein_count = len(prots)
		for p in prots:
			self.min_genome_length += len(self.proteome_seq[p])
			if p not in self.protein_to_accs:
				discard = self.proteome_seq.pop(p)
		#aa to nt
		self.min_genome_length = self.min_genome_length * 3
	
	def convert_prots_to_tetras(self):
		for p in self.proteome_seq:
			tetras, seqlen, valid_kmer_ct = self.tt.tetras_to_base20(self.proteome_seq[p])
			self.proteome_seq[p] = tetras
			self.kmer_cts[p] = valid_kmer_ct
			self.seqlens[p] = seqlen
			
		
	def jsonify(self):
		json_records = {}
		json_records['genome'] = self.basename
		json_records['genome_length'] = self.genome_length
		json_records['coding_bases'] = self.min_genome_length
		json_records['protein_count'] = self.protein_count
		json_records['taxonomy'] = self.taxonomy
		json_records['SCPs'] = {}
		
		if len(self.protein_to_accs) >= 5:
		
			for p in self.protein_to_accs:
				json_records['SCPs'][p] = {"hmm_accession":self.protein_to_accs[p],
									"hmm_score":self.accs_to_besthits[self.protein_to_accs[p]],
									"protein_len":self.seqlens[p],
									"valid_tetramer_ct":self.kmer_cts[p],
									"base_20_tetras":self.proteome_seq[p].tolist()}
									
			json_records = json.dumps(json_records, indent = 4)
			self.json = json_records
		else:
			self.json_records = None
			
	def set_values_from_mem(self, basename, genlength = 0, p_to_acc = None, acc_to_bh = None, prot_seqdict = None):
		self.basename = basename
		self.genome_length = genlength
		self.protein_to_accs = p_to_acc
		self.accs_to_besthits = acc_to_bh
		self.proteome_seq = prot_seqdict
	
	def set_values_from_files(self, protein_file, hmm_file, genome_file = None):
		self.pf = protein_file
		self.hf = hmm_file
		self.gf = genome_file
		self.basename = get_basename(self.gf)
		self.load_hmm()
		if self.protein_to_accs is not None:
			self.load_genome()
			self.load_prots()

	def jsonify_loaded_data(self):
		self.json = None
		if self.protein_to_accs is not None and self.proteome_seq is not None:
			self.filter_prots()
			self.convert_prots_to_tetras()
			self.jsonify()
			
			#Clean up
			self.protein_to_accs = None
			self.accs_to_besthits = None
			self.proteome_seq = None
			self.seqlens = {}
			self.kmer_cts = {}		
		
		return self.json
		
	def write_crystal(self, output):
		if self.compress:
			content = self.json.encode(encoding = "ascii")
			with gzip.open(output+".gz", "wb") as outwriter:
				outwriter.write(content)
			
		else:
			with open(output, "w") as outwriter:
				outwriter.write(self.json)
	
	
'''
import os
global d1
global d2
global d3

genome_files = None

protdir = sys.argv[1]
hmmdir = sys.argv[2]


if len(sys.argv) > 3:
	genomedir = sys.argv[3]
	genome_files = os.listdir(genomedir)
	genome_files.sort()
	genome_files = [genomedir + f for f in genome_files]

hmm_files = os.listdir(hmmdir)
protein_files = os.listdir(protdir)

hmm_files.sort()
protein_files.sort()

protein_files = [protdir + f for f in protein_files]
hmm_files = [hmmdir + f for f in hmm_files]

d1 = {}
d2 = {}
d3 = None

for p in protein_files:
	pbase = get_basename(p)
	#pbase = pbase.split(".aa_genes.txt")[0]
	pbase = pbase.split(".genes.fasta")[0]
	
	d1[pbase] = p
	#print(pbase)
	
for h in hmm_files:
	hbase = get_basename(h)

	hbase = hbase.split(".hmmsearch.txt")[0]
	#hbase = hbase.split("_SCP_HMM_search")[0]
	hbase = hbase.replace("fastaai_v2_prokaryotes_", "")
	
	d2[hbase] = h
	
if genome_files is not None:
	d3 = {}
	for g in genome_files:
		gbase = get_basename(g)
		if gbase not in d1:
			gbase = gbase.replace("_genomic", "")
			if gbase not in d1:
				gbase = gbase.split(".fasta")[0]
			
		d3[gbase] = g

def run_prep(p):
	mn = crystalizer()	
	if p in d1 and p in d2:
		pfile = d1[p]
		hfile = d2[p]
		
		if d3 is not None:
			if p in d3:
				gfile = d3[p]
			else:
				print("Genome for protein", p, "missing. It will be skipped.")
				gfile = None
		else:
			gfile = None
			
		j = mn.prep_hmm(pfile, hfile, gfile)
		
		if j is not None:
			with open("filtered_crystals/"+p+"_cyrstal.txt", "w") as fh:
				fh.write(j)
		else:
			print("Genome", p, "did not recover at least 5 single-copy proteins. It cannot be used by FastAAI.")
		
	
	return p


if not os.path.exists("filtered_crystals/"):
	os.mkdir("filtered_crystals/")

protnames = list(d1.keys())
#protnames = protnames[0:100]

ct = 0
pool = multiprocessing.Pool(20)
for res in pool.imap_unordered(run_prep,protnames):
	#print(ct)
	ct += 1

pool.close()
pool.join()
'''




	