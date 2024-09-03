import sys
import os


from fastaai_db_prep import initialize_fastaai_database
from fastaai_preprocessor import fastaai_preprocess
from fastaai_consumer import fastaai_consume
from fastaai_final_db_cleanup import finalize_db
from fastaai_query import run_fastaai_query

def main():
	valid_options = ['init',
					'preproc',
					'consume',
					'finalize',
					'query']
	
	if len(sys.argv) < 2:
		print("Please select a module from one of the following:")
		for m in valid_options:
			print("\t", m)
		sys.exit()
		
	module = sys.argv[1]
	if module not in valid_options:
		print("Module not recognized")
		print("Please select a module from one of the following:")
		for m in valid_options:
			print("\t", m)
			
		sys.exit()
		
	
		
	if module == 'init':
		initialize_fastaai_database()
	if module == 'preproc':
		fastaai_preprocess()
	if module == 'consume':
		fastaai_consume()
	if module == 'finalize':
		finalize_db()
	if module == 'query':
		run_fastaai_query()
		

main()