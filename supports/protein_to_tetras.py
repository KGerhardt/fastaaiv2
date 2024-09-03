import numpy as np

class tetra_mapper:
	def __init__(self):
		self.aa_to_idx = {"A":0,
						"C":1,
						"D":2,
						"E":3,
						"F":4,
						"G":5,
						"H":6,
						"I":7,
						"K":8,
						"L":9,
						"M":10,
						"N":11,
						"P":12,
						"Q":13,
						"R":14,
						"S":15,
						"T":16,
						"V":17,
						"W":18,
						"Y":19}

	def tetras_to_base20(self, seq):
		seqlen = len(seq)
		#num tetramers = len(seq) - 4 + 1, just make it -3.
		n_kmers = seqlen - 3
		
		#Converts the characters in a sequence into their ascii int value
		as_ints = []
		for i in seq:
			if i in self.aa_to_idx:
				next_char = self.aa_to_idx[i]
			else:
				next_char = -1
			
			as_ints.append(next_char)
		#as_ints = np.array([self.aa_to_idx(i) for i in seq if i in self.aa_to_idx], dtype = np.int32)
		as_ints = np.array(as_ints, dtype = np.int32)
		
		#create seq like 0,1,2,3; 1,2,3,4; 2,3,4,5... for each tetramer that needs a value
		kmers = np.arange(4*n_kmers)
		kmers = kmers % 4 + kmers // 4
		
		#Select the characters (as ints) corresponding to each tetramer all at once and reshape into rows of 4, 
		#each row corresp. to a successive tetramer
		kmers = as_ints[kmers].reshape((n_kmers, 4))
		
		#Remove illegal characters not in the aa_to_idx dict
		kmers = kmers[np.all(kmers >= 0, axis = 1)]
		
		
		#Given four 2-digit numbers, these multipliers work as offsets so that all digits are preserved in order when summed
		mult = np.array([8000, 400, 20, 1], dtype = np.int32)

		#Convert tetramers to base 20.
		kmers =  np.dot(kmers, mult)
		
		#unique, but preserve ordering
		_, idx = np.unique(kmers, return_index=True)
		kmers = kmers[np.sort(idx)]
		
		#Has to come after unique or it screws up jacc. calculations later.
		num_kmers = kmers.shape[0]
		
		return kmers, seqlen, num_kmers
	


