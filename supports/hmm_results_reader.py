try:
	from .agnostic_reader import agnostic_reader
except:
	from agnostic_reader import agnostic_reader

import numpy as np

class hmm_results:
	def __init__(self, results_file):
		self.f = results_file
		
		self.contents = None
		
		self.protein_to_acc = None
		self.acc_to_score = None
		
	def load_and_besthit(self, best_per_scp = True, reciprocal_besthit = True, filter_to = None):
		hmm_file = np.loadtxt(self.f, comments = '#', delimiter = "\t", dtype = str, ndmin=2)
		if hmm_file.shape[1] > 5: #Doesn't matter how many cols as long as the number is >1
			#Columns
			#0 = protein
			#2 = maybe acc?
			#3 = accession
			#8 = score
			#hmm_file = np.loadtxt(hmm_file_name, comments = '#', usecols = (0, 3, 8), dtype=(str))
			#Sort the hmm file based on the score column in descending order.
			#hmm_file = hmm_file[hmm_file[:,2].astype(float).argsort()[::-1]]
			#hmm_file = hmm_file[hmm_file[:,8].astype(float) >= 30.0] #busco
			
			#Filter to some set of SCPs before best hitting
			if filter_to is not None:
				ok_rows = []
				row = 0
				for scp_acc in hmm_file[:,3]:
					if scp_acc in filter_to:
						ok_rows.append(row)
					row += 1
				ok_rows = np.array(ok_rows, dtype = np.int32)
				hmm_file = hmm_file[ok_rows]				

			hmm_file = hmm_file[hmm_file[:,8].astype(float).argsort()[::-1]]
			
			if best_per_scp:
				#Filter the file again for the unique ACCESSION names, since we're only allowed one gene per accession, I guess?
				#Don't sort the indices, we don't care about the scores anymore.
				hmm_file = hmm_file[np.unique(hmm_file[:,3], return_index = True)[1]]

			if reciprocal_besthit:
				#Identify the first row where each gene name appears, after sorting by score; 
				#in effect, return the highest scoring assignment per gene name
				#Sort the indices of the result to match the score-sorted table instead of alphabetical order of gene names
				
				#Filter per gene
				hmm_file = hmm_file[np.sort(np.unique(hmm_file[:,0], return_index = True)[1])]
				
		else:
			hmm_file = None
		
		if hmm_file is not None:
			if len(hmm_file) > 0:
				self.contents = hmm_file
		
		
	def format(self):
		if self.contents is not None:
			self.protein_to_acc = {}
			self.acc_to_score = {}
			idx = 0
			for prot in self.contents[:,0]:
				acc = self.contents[idx, 3]
				#acc = self.contents[idx, 2] #busco
				self.protein_to_acc[prot] = acc
				score = round(float(self.contents[idx, 8]), 6)
				self.acc_to_score[acc] = score
				idx += 1
				
				
	def run(self, best_per_scp = True, reciprocal_besthit = True, filter_to = None):
		self.load_and_besthit(best_per_scp = best_per_scp, reciprocal_besthit = reciprocal_besthit, filter_to = filter_to)
		self.format()
				
		return self.protein_to_acc, self.acc_to_score
				
	
		
			