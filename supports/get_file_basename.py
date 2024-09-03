import os

def get_basename(filepath):
	fn = os.path.basename(filepath)
	#Special gzip case
	if fn.endswith(".gz"):
		fn = os.path.splitext(fn)[0]
	
	fn = os.path.splitext(fn)[0]
	
	if '.genes.fasta' in fn:
		fn = fn.split(".genes.fasta")[0]
	if '.hmmsearch.txt' in fn:
		fn = fn.split(".hmmsearch.txt")[0]
	
	return fn