import gzip
import os

from .agnostic_reader import agnostic_reader

#Plain + gzip fasta file importer
#Returns dict of seqid:{'description': , 'sequence':}
class fasta_file:
	def __init__(self, file):
		self.file_path = file
		self.name = None
		#Guarantees same order as the file was read
		self.sequence_name_order = []
		#Contains the file's seq ids, descriptions, and sequences
		#self.contents = []
		self.deflines = {}
		self.sequences = {}
		
	def get_file_basename(self):
		filename = os.path.basename(self.file_path)
		self.name = filename
		
	def file_contents_to_tuples(self, fasta_file_contents):
		ordered_results = []
				
		#Split on seqid starts; removes these characters from the result set
		#first item of the resulting list is always empty string '', so remove it
		#Returns a list of contigs and their sequences
		fasta_file_contents = fasta_file_contents.split(">")[1:]
		
		#For each contig:
		for sequence in fasta_file_contents:
			#The first line is always the defline; first newline char indicates its end
			sequence = sequence.split("\n", maxsplit=1)
			#Deflines consist of >seqid[whitespace]description;
			#'>' has already been removed, so split on first whitespace and return halves
			seqline = sequence[0].split(maxsplit = 1)
			seqid = seqline[0]
			#self.sequence_name_order.append(seqid)
			
			if len(seqline) > 1:
				desc = seqline[1]
			else:
				desc = ""
				
			#The remaining sequence is all the lines after the first, 
			#so remove newlines
			sequence = sequence[1].replace('\n', '')
			
			next_row = (seqid, desc, sequence,)
			ordered_results.append(next_row)
			
		return ordered_results
		
	def import_fasta(self):
		self.contents = []
		self.sequence_name_order = []
		
		self.get_file_basename()
		
		ar = agnostic_reader(self.file_path)
		fasta_file_contents = ar.read()

		fasta_file_contents = self.file_contents_to_tuples(fasta_file_contents)
				
		#For each contig:
		for fasta_tuple in fasta_file_contents:
			seqid = fasta_tuple[0]
			desc = fasta_tuple[1]
			sequence = fasta_tuple[2]
			self.sequence_name_order.append(seqid)
			self.deflines[seqid] = desc
			self.sequences[seqid] = sequence
			
		return None
		
	#Easy function to extract the file to another program.
	def read(self):
		self.sequence_name_order = []
		self.deflines = {}
		self.sequences = {}
		self.import_fasta()
		
	def tuples_to_file(self):
		fasta_file = None
		if len(self.sequence_name_order) > 0:
			fasta_file = []
			for seqid in self.sequence_name_order:
				desc = self.deflines[seqid]
				sequence = self.sequences[seqid]
				next_line = seqid + " " + desc + "\n" + sequence
				fasta_file.append(next_line)
				
			fasta_file = '\n'.join(fasta_file)
			
		self.contents = fasta_file