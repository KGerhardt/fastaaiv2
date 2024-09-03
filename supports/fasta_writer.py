import gzip

class fasta_writer:
	def __init__(self):
		self.contents = None
	
	#Pyrodigal standard: nt = 70 chars/line max, aa = 60
	def num_bp_line_format(self, string, size = 60):
		#ceiling funciton without the math module
		ceiling = int(round((len(string)/size)+0.5, 0))
		formatted = '\n'.join([string[(i*size):(i+1)*size] for i in range(0, ceiling)])
		#remove occasional trailing newline
		formatted = formatted.strip()
		return formatted
	
	def info_to_fasta(self, results_tuples, protein = False):
		self.contents = []
		for tup in results_tuples:
			defline = ">" + tup[0] + " " + tup[1]
			if protein:
				chars_per_line = 60
			else:
				chars_per_line = 70
				
			formatted_seq = self.num_bp_line_format(tup[2], size = chars_per_line)
				
			self.contents.append(defline)
			self.contents.append(formatted_seq)
		
		#We need an extra char so the file actually ends with a newline
		self.contents.append("")
		self.contents = '\n'.join(self.contents)
		
			
	#Also functions as a general gzip-writer
	def write_file(self, output_path = None, is_gz = False, compression_level = 6):
		if output_path is not None and self.contents is not None:
			if is_gz:
				self.contents = self.contents.encode(encoding = 'utf-8')
				self.contents = gzip.compress(self.contents, 
											compresslevel = compression_level,
											mtime=0)
				if not output_path.endswith(".gz"):
					output_path += '.gz'
					
				with open(output_path, mode = 'wb') as fh:
					fh.write(self.contents)
					
			else:
				with open(output_path, 'w') as fh:
					fh.write(self.contents)
			
		#Reset
		self.contents = None
		
		