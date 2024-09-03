import sys
import os
import argparse
import multiprocessing
import sqlite3

from supports.agnostic_reader import agnostic_reader
from supports.fasta_importer import fasta_file

from supports.get_file_basename import get_basename

from supports.pyrodigal_manager import new_pyrodigal_manager
from supports.pyhmmer_manager import new_pyhmmer_manager
from supports.crystalizer import crystalizer
from supports.fasta_importer import fasta_file

class fastaai_preprocessor:
	def __init__(self, dbpath, genome_directory, module = "prokaryotes", output_nt = False, threads = 1, compress = True):
		self.dbpath = dbpath
		self.input_genomes = genome_directory
		self.module = module

		self.tax_by_genome = {}

		self.output_nt = output_nt

		self.threads = threads
		self.compress = compress
		
	def make_directories(self):
		if not os.path.exists("proteins_aa"):
			os.mkdir("proteins_aa")
		if not os.path.exists("proteins_nt"):
			os.mkdir("proteins_nt")
		if not os.path.exists("hmm_searches"):
			os.mkdir("hmm_searches")
		if not os.path.exists("crystals"):
			os.mkdir("crystals")
		
	def load_taxonomy(self, taxfile):
		self.tax_by_genome = {}
		with open(taxfile) as fh:
			for line in fh:
				segs = line.strip().split("\t")
				genome = segs[0]
				gtdb_fmt_tax_string = segs[1]
				self.tax_by_genome[genome] = gtdb_fmt_tax_string
		
	def run_preprocessing(self):
		genome_files = os.listdir(self.input_genomes)
		genome_files = [os.path.normpath(self.input_genomes + "/" + g) for g in genome_files]
		if self.threads > len(genome_files):
			self.threads = len(genome_files)
			
		for g in genome_files:
			bn = get_basename(g)
			if bn not in self.tax_by_genome:
				self.tax_by_genome[bn] = 'd__;p__;c__;o__;f__;g__;s__'
			
		self.make_directories()
		
		if self.module == "prokaryotes":
			pyhmmerer = new_pyhmmer_manager(self.compress)
			pyhmmerer.load_hmm_from_database(self.dbpath)
			
			preproc_args = [prokaryotic_input_file(input_genome = g,
													taxonomy = self.tax_by_genome[get_basename(g)],
													compress = self.compress,
													pyhmmer_manager = pyhmmerer) for g in genome_files]
			
		if self.module == "eukaryotes":
			pass
			#Still working on this. Need to decide on metaeuk or AUGUSTUS. AUGUSTUS sucks to install.

		pool = multiprocessing.Pool(self.threads)
		for result in pool.imap_unordered(run_preprocessor, preproc_args):
			print(result)
			
		pool.close()
		pool.join()


class prokaryotic_input_file:
	def __init__(self, 
				genome_type = "prokaryote", 
				preprocessing_stage = "genome",
				add_rules = None,
				input_genome = None,
				input_protein_aa = None,
				input_protein_nt = None,
				input_hmm = None,
				input_crystal = None,
				taxonomy = 'd__;p__;c__;o__;f__;g__;s__',
				compress = False,
				nucleotide_output = False,
				prot_nt_out = "proteins_nt",
				prot_aa_out = "proteins_aa",
				hmm_out = "hmm_searches",
				crystal_out = "crystals",
				pyhmmer_manager = None
				):
		
		self.gtype = genome_type
		self.start_point = preprocessing_stage
		self.genome_class_rules = add_rules
		
		self.genome_file = input_genome
		self.protein_file_aa = input_protein_aa
		self.protein_file_nt = input_protein_nt
		
		self.hmm_file = input_hmm
		self.crystal_file = input_crystal
		
		self.name_base = get_basename(self.genome_file)
		
		self.genome = None
		self.protein = None
		self.protein_nt_ct = None
		self.hmm_best_hits = None
		self.hmm_acc_to_score = None
		self.crystal = None
		
		self.prot_aa_out_dir = prot_aa_out
		self.prot_nt_out_dir = prot_nt_out
		self.hmm_out_dir = hmm_out
		self.crystal_out_dir = crystal_out
		
		self.compress = compress
		self.output_nt = nucleotide_output
		
		self.pyhmmer_manager = pyhmmer_manager
		
		self.genome_length = 0
		self.num_proteins = 0
		self.num_scps = 0
		
		self.gtdb_tax_string = taxonomy
		
	def genome_to_protein(self):
		if self.genome_file is None and self.protein_file is None:
			print("FastAAI cannot proceed without either a genome or a protein file")
			return None
	
		if self.genome_file is not None:
			fr = fasta_file(self.genome_file)
			fr.read()
			for seq in fr.sequences:
				self.genome_length += len(fr.sequences[seq])
				
			self.genome = fr.sequences
			fr = None
		
		if self.protein_file_aa is None:
			self.protein_file_aa = os.path.normpath(self.prot_aa_out_dir + "/" +self.name_base+"_protein_aa.txt")
			if self.output_nt:
				self.protein_file_nt = os.path.normpath(self.prot_aa_out_dir + "/" +self.name_base+"_protein_nt.txt")
			else:
				self.protein_file_nt = None
			
		
		mn = new_pyrodigal_manager()
		self.protein, aa_filepath_copy = mn.run_in_mem(genome_dict = self.genome, compress = self.compress, outnt = self.protein_file_nt, outaa = self.protein_file_aa)
		self.protein_nt_ct = mn.protein_length_nt
		mn = None	
			
	def protein_to_hmm(self):
		if self.protein is None:
			if self.protein_file is None:
				print("Cannot load a nonexistent protein file.")
			else:
				fr = fasta_file(self.protein_file)
				fr.read()
				self.protein = fr.seqeunces
				fr = None
				
		if self.protein is not None:
			self.pyhmmer_manager.convert_protein_seqs_in_mem(self.protein)
			self.pyhmmer_manager.search_protein()
			self.hmm = os.path.normpath(self.hmm_out_dir + "/" + self.name_base + "_hmmsearch.txt")
			self.pyhmmer_manager.to_hmm_file_pyhmmer_method(self.hmm)
			self.hmm_best_hits, self.hmm_acc_to_score = self.pyhmmer_manager.filter_to_best_hits()
		
	def run_crystalize(self):
		mn = crystalizer(taxonomy = self.gtdb_tax_string, compress = self.compress)
		mn.set_values_from_mem(basename = self.name_base,
								genlength = self.genome_length, 
								p_to_acc = self.hmm_best_hits, 
								acc_to_bh = self.hmm_acc_to_score, 
								prot_seqdict = self.protein)
								
		self.crystal = mn.jsonify_loaded_data()
		if self.crystal_file is None:
			self.crystal_file = os.path.normpath(self.crystal_out_dir + "/" +self.name_base+"_fastaai_crystal.txt")
		
		mn.write_crystal(self.crystal_file)
		
	def run(self):
		self.genome_to_protein()
		self.protein_to_hmm()
		self.run_crystalize()
		return self.name_base
	
	
def hmm_preproc_initializer(hmm_database, do_compress = False):
	global pyhmmer_manager
	pyhmmer_manager = new_pyhmmer_manager(do_compress)
	pyhmmer_manager.load_hmm_from_database(hmm_database)
	return pyhmmer_manager
	
def run_preprocessor(preproc_obj):
	this_genome = preproc_obj.run()
	return this_genome
	
	
import argparse
def consumer_options():
	parser = argparse.ArgumentParser(description='Process genomes into FastAAI v2 crystals. Needs a database that has been initialized with fastaai2 init and a set of genomes')
	
	parser.add_argument('-db', '--database', dest='database', default = None,
						help='An initialized FastAAI v2 database')
	parser.add_argument('-g', '--genomes', dest='genomes_dir', default = None,
						help='A directory containing genomes to preprocess')
	parser.add_argument('--taxonomy_file', dest='tax', default = None,
						help='A tab-separated file containing genome ID in column 1, GTDB-formatted taxonomy string in column 2. No header.')
	parser.add_argument('--eukaryotic', dest='is_euk', action = 'store_true',
						help='Flag to indicate that genomes are eukaryotic.')						
	parser.add_argument('--output_nt', dest='nt', action = 'store_true',
						help='Flag for producing nucleotide sequences for predicted proteins in addition to amino acid sequences. These are not used by FastAAI, this is just for your convenience.')						
	parser.add_argument('--compress', dest='compress', action = 'store_true',
						help='Flag for gzipping outputs. Modestly increases runtime.')						
											
	parser.add_argument('-t', '--threads', dest = 'threads', default = 1, type = int,
						help='How many processes to use')


	args = parser.parse_known_args()
	return parser, vars(args[0])

def fastaai_preprocess():
	parser, args = consumer_options()
	if len(sys.argv) < 3:
		parser.print_help()
	else:
		db = args['database']
		genomes = args['genome_dir']
		tax_file = args['tax']
		euk = args['is_euk']
		do_nt = args['nt']
		comp = args['compress']
		threads = args['threads']
		
		if euk:
			org = 'eukaryotes'
		else:
			org = 'prokaryotes'

		if not os.path.exists(db):
			print("A FastAAI v2 database must already exist at path", db, "or this step cannot proceed.")
			print("Please check that a database already exists or make one with fastaai2 init")
		else:
			mn = fastaai_preprocessor(dbpath = db, genome_directory = genomes, module = org, output_nt = do_nt, threads = threads, compress = comp)
			mn.load_taxonomy(tax_file)
			mn.run_preprocessing()	