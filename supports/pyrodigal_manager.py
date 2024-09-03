import sys
import os

import pyrodigal

from .fasta_importer import fasta_file
import io
import gzip

from .get_file_basename import get_basename

import multiprocessing

class new_pyrodigal_manager:
	def __init__(self, trans_tables = [11, 4], meta = False, verbose = False):
		self.genome_dict = None
		self.training_seq = None
		
		self.genome_length = 0
		
		self.meta = meta
		self.verbose = verbose
		
		self.training_seq = None
		self.gene_finder = None
		
		self.coding_bases_by_table = {}
		self.coding_density_by_table = {}
		
		self.trans_tables = trans_tables
		
		self.current_genes = None
		self.winning_coding_dens = 0
		self.winning_table = None
		
		self.genes_nt = None
		self.genes_aa = None
		
		self.protein_seqs = {}
		self.protien_length_nt = 0
		
	#Reset to initial conditions so that the current use of this object doesn't interfere with the next
	def reset(self):
		self.genome_dict = None
		self.training_seq = None
		
		self.genome_length = 0
		
		self.training_seq = None
		self.gene_finder = None
		
		self.coding_bases_by_table = {}
		self.coding_density_by_table = {}
				
		self.current_genes = None
		self.winning_coding_dens = 0
		self.winning_table = None
		
		self.genes_nt = None
		self.genes_aa = None
		
		self.protein_seqs = {}
		self.protein_length_nt = 0
		
		
	def load_genome_from_file(self, genome_file):
		ff = fasta_file(genome_file)
		ff.import_fasta()
		
		self.prep_genome_dict_for_prediction(ff.sequences)
		
	def prep_genome_dict_for_prediction(self, genome_dict):
		self.reset()
		
		self.genome_dict = {}
		for g in genome_dict:
			#asbytes = g.encode(encoding = "ascii")
			self.genome_length += len(genome_dict[g])
			self.genome_dict[g] = genome_dict[g].encode(encoding = "ascii")
		
	def prep_training_seq(self):
		running_sum = 0
		seqs_added = 0
		if self.training_seq is None:
			self.training_seq = []
			for seq in self.genome_dict:
				running_sum += len(self.genome_dict[seq])
				if seqs_added > 0:
					#Prodigal interleaving logic - add this breaker between sequences, starting at sequence 2
					self.training_seq.append(b'TTAATTAATTAA')
					running_sum += 12
					
				seqs_added += 1
					
				#Handle excessive size
				if running_sum >= 32000000:					
					print("Warning:  Sequence is long (max 32000000 for training).")
					print("Training on the first 32000000 bases.")
				
					to_remove = running_sum - 32000000
					
					#Remove excess characters
					cut_seq = self.genome_dict[seq][:-to_remove]
					#Add the partial seq
					self.training_seq.append(cut_seq)
					
					#Stop the loop and move to training
					break
				
				#add in a full sequence
				self.training_seq.append(self.genome_dict[seq])

			if seqs_added > 1:
				self.training_seq.append(b'TTAATTAATTAA')
				
			self.training_seq = b''.join(self.training_seq)
		
		if len(self.training_seq) < 20000:
			if self.verbose:
				print("Can't train on 20 thousand or fewer characters. Switching to meta mode.")
			self.manager = pd.GeneFinder(meta=True)
			self.meta = Trues
		else:
			if self.verbose:
				#G is 71, C is 67; we're counting G + C and dividing by the total.
				gc = round(((self.training_seq.count(67) + self.training_seq.count(71))/ len(self.training_seq)) * 100, 2)
				print(len(self.training_seq), "bp training seq created,", gc, "pct GC")
				
	def train(self, table = 11):
		if not self.meta:
			self.gene_finder = pyrodigal.GeneFinder(meta = False)
			self.gene_finder.train(self.training_seq, translation_table = table)
		else:
			self.gene_finder = pyrodigal.GeneFinder(meta = True)
			
	def predict(self):
		current_table = self.gene_finder.training_info.translation_table
		
		if current_table is None:
			current_table = "meta"
			
		self.coding_bases_by_table[current_table] = 0
		
		next_gene_set = {}
		
		if self.gene_finder is not None:
			for seq in self.genome_dict:
				genes = self.gene_finder.find_genes(self.genome_dict[seq])
				next_gene_set[seq] = genes
				for g in genes:	
					self.coding_bases_by_table[current_table] += len(g.sequence())
		
		self.coding_density_by_table[current_table] = self.coding_bases_by_table[current_table] / self.genome_length
		
		
		if self.current_genes is None:
			self.current_genes = next_gene_set
			self.winning_coding_dens = self.coding_density_by_table[current_table]
			self.winning_table = current_table
		else:
			if self.coding_density_by_table[current_table] > (self.winning_coding_dens * 1.1):
				if self.verbose:
					print("Translation table", self.winning_table, "had coding density of", self.winning_coding_dens)
					print("New translation table", current_table, "more than 10% better, coding density of", self.coding_density_by_table[current_table])
					print("New table's predicted genes will be used")
				self.current_genes = next_gene_set
				self.winning_coding_dens = self.coding_density_by_table[current_table]
				self.winning_table = current_table
			
		
	#Format as fastaai expects / as would be loaded by supports/fasta_importer
	def genes_to_mem(self):
		self.protein_seqs = {}
		for seq in self.current_genes:
			seqid = 1
			for g in self.current_genes[seq]:
				gene_length = (g.end - g.begin + 1)
				self.protein_length_nt += gene_length
				gene_name = seq + "_" + str(seqid)
				seqid += 1
				self.protein_seqs[gene_name] = g.translate(None)
			
	def write_nt(self, compress = False, outfile = None):
		if outfile is not None:
			self.genes_nt = io.StringIO("")
			sz = 0
			for seq in self.current_genes:
				sz += len(self.current_genes[seq])
				self.current_genes[seq].write_genes(self.genes_nt, seq)
				
			self.genes_nt = self.genes_nt.getvalue()
			
			if compress:
				self.genes_nt = self.genes_nt.encode(encoding = "ascii")
				self.genes_nt = gzip.compress(self.genes_nt)
				with open(outfile+".gz", 'wb') as out:
					out.write(self.genes_nt)
			else:
				with open(outfile, "w") as out:
					out.write(self.genes_nt)
						
			#This never needs returned in fastaai
			self.genes_nt = None
		
		
	def write_aa(self, compress = False, outfile = None):
		if outfile is not None:
			self.genes_aa = io.StringIO("")
			
			for seq in self.current_genes:
				self.current_genes[seq].write_translations(self.genes_aa, seq)
				
			self.genes_aa = self.genes_aa.getvalue()
		
			if compress:
				self.genes_aa = self.genes_aa.encode(encoding = "ascii")
				with gzip.open(outfile+".gz", 'wb') as out:
					out.write(self.genes_aa)
			else:
				with open(outfile, "w") as out:
					out.write(self.genes_aa)
					
			self.genes_aa = None

	
	def run(self, genome_file, compress = False, outnt = None, outaa = None):
		self.load_genome_from_file(genome_file)
		self.prep_training_seq()
		for table in self.trans_tables:
			self.train(table = table)
			self.predict()
			
		self.genes_to_mem()
		
		self.write_nt(compress = compress, outfile = outnt)
		self.write_aa(compress = compress, outfile = outaa)
				
		return self.protein_seqs, outaa
		
	def run_in_mem(self, genome_dict, compress = False, outnt = None, outaa = None):
		self.prep_genome_dict_for_prediction(genome_dict)
		self.prep_training_seq()
		for table in self.trans_tables:
			self.train(table = table)
			self.predict()
			#print(len(mn.current_genes))
			
		self.genes_to_mem()
		
		self.write_nt(compress = compress, outfile = outnt)
		self.write_aa(compress = compress, outfile = outaa)
				
		return self.protein_seqs, outaa
		
	def run_for_fastaai():
		pass

def run_job(args):
	genome = args[0]
	output_nt = args[1]
	output_aa = args[2]
	return_fp = args[3]
	compress = args[4]
	compared_tables = args[5]
	
	mn = new_pyrodigal_manager(trans_tables = compared_tables)
	protein_seq_dict = mn.run(genome_file = genome, compress = compress, outnt = output_nt, outaa = output_aa)
	
	if return_fp:
		protein_seq_dict = output_aa
	
	return genome, protein_seq_dict


def operate_in_parallel(genome_files, output_nt = None, output_aa = None, return_file_path_instead_of_protein_dict = False, 
						threads = 1, compress = True, trans_tabs = [11, 4]):
	
	if output_nt is None:
		output_nt = [None] * len(genome_files)
	if output_aa is None:
		output_aa = [None] * len(genome_files)
		
	argset = []
	for g, n, a in zip(genome_files, output_nt, output_aa):
		next_arg = (g, n, a, return_file_path_instead_of_protein_dict, compress, trans_tabs,)
		argset.append(next_arg)
	
	output_set = {}
	pool = multiprocessing.Pool(threads)
	for result in pool.imap_unordered(run_job, argset):
		og_file = result[0]
		seqdict_or_filepath = result[1]
		output_set[og_file] = seqdict_or_filepath
		
	pool.close()
	pool.join()
		
'''		
genome_dir = sys.argv[1]
threads = int(sys.argv[2])
		
genomes = os.listdir(genome_dir)
genomes = [os.path.normpath(genome_dir+"/"+g) for g in genomes]
aas = []
nts = []

if not os.path.exists("prot_aa"):
	os.mkdir("prot_aa")
	
if not os.path.exists("prot_nt"):
	os.mkdir("prot_nt")

mn = new_pyrodigal_manager(trans_tables = [11, 4])
for g in genomes:
	prf = get_basename(g)
	aa = "prot_aa/"+prf+".aa_genes.txt"
	aas.append(aa)
	nt = "prot_nt/"+prf+".nt_genes.txt"
	nts.append(nt)


operate_in_parallel(genome_files = genomes, output_aa = aas, threads = 20)
'''

