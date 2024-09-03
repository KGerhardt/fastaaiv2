import sys
import os
try:
	from .agnostic_reader import agnostic_reader
except:
	from agnostic_reader import agnostic_reader

try:
	from .get_file_basename import get_basename
except:
	from get_file_basename import get_basename

#Augustus block profiles happen one at a time instead of a whole SCP search, I think
def read_block_hmm(file):
	file_contents = []
	block_names = []
	file_name = get_basename(file)
	file_name = file_name.replace(".hmm", "")
	scp_name = None
	take_next = False
	fh = agnostic_reader(file)
	contents = fh.read()
	contents = contents.splitlines()
	for line in contents:
		file_contents.append(line)
		if take_next:
			scp_name = line.strip()
			take_next = False
		if line.startswith("[name]"):
			take_next = True
		if line.startswith("name="):
			next_block_name = line.strip()[5:]
			block_names.append(next_block_name)
		
	file_contents = '\n'.join(file_contents)
	
	return scp_name, file_name, file_contents
	
def read_euk_hmm_directory(d):
	hmm_model_by_acc, acc_to_name = {}, {}
	for f in os.listdir(d):
		scp, file, cont = read_block_hmm(os.path.normpath(d+"/"+f))
		hmm_model_by_acc[file] = cont
		acc_to_name[file] = file
		
	return hmm_model_by_acc, acc_to_name