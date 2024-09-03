import sys
import os
import numpy as np
import sqlite3 
import json

from supports.hmm_model_reader import read_prok_hmm_directory
from supports.augustus_block_hmm_reader import read_euk_hmm_directory
from supports.hmm_jsonifier import hmm_jsonifier

class fastaai_db_prepper:
	def __init__(self, path):
		self.path = path
		self.conn = None
		self.curs = None
		
	
	def open(self):
		if self.conn is None:
			self.conn = sqlite3.connect(self.path)
			self.curs = self.conn.cursor()
	
	def close(self):
		if self.conn is not None:
			self.curs.close()
			self.conn.close()
			self.conn = None
			self.curs = None
		
	def initialize_genome_index_table(self):
		
		
		gclass_fields = ', '.join([
				"genome_class INTEGER PRIMARY KEY",
				"class_label TEXT",
				'scp_members',
				"rationale TEXT"])
		
	
		metadata_fields = ', '.join([
				"genome_name TEXT PRIMARY KEY",
				"genome_id INTEGER",
				"genome_length INTEGER",
				"proteome_ct INTEGER",
				'scp_count INTEGER',
				'genome_class INTEGER',
				'genome_class_label TEXT',
				'dom TEXT',
				'phyl TEXT',
				'cl TEXT',
				'ord TEXT',
				'fam TEXT',
				'genu TEXT',
				'spec TEXT'
				])

		
		sql = "CREATE TABLE IF NOT EXISTS genome_classes ({fields})".format(fields = gclass_fields)
		self.curs.execute(sql)
		sql = "INSERT OR IGNORE INTO genome_classes VALUES (?, ?, ?, ?)"
		self.curs.execute(sql, (0, "Default genome class", None, "Any otherwise unlabeled genome is identified as being part of class zero",))
		
		
		sql = "CREATE TABLE IF NOT EXISTS genome_metadata ({fields})".format(fields = metadata_fields)
		self.curs.execute(sql)
		
	
		sql = "CREATE INDEX IF NOT EXISTS gclass_idx ON genome_classes (genome_class)"
		self.curs.execute(sql)
		
		
		sql = "CREATE INDEX IF NOT EXISTS genomes_index ON genome_metadata (genome_id)"
		self.curs.execute(sql)
		
		self.conn.commit()
			
	def initialize_gak(self, scp_dir, gen_type = "prokaryote"):
		scp_fields = ', '.join([
				"scp_name TEXT",
				"scp_acc TEXT PRIMARY KEY",
				"organism_type TEXT",
				"scp_id INTEGER",
				"scp_model TEXT"
				])
	
		gak_fields = ', '.join([
				"genome_id INTEGER",
				"scp_id INTEGER",
				"protein_length INTEGER",
				"hmm_score REAL",
				"tetra_count INTEGER",
				"PRIMARY KEY (genome_id, scp_id)"
				])
		
		sql = "CREATE TABLE IF NOT EXISTS scp_data ({fields})".format(fields = scp_fields)
		self.curs.execute(sql)
		
		sql = "CREATE TABLE IF NOT EXISTS genome_accession_counts ({fields})".format(fields = gak_fields)
		self.curs.execute(sql)
		
		
		sql = "CREATE INDEX IF NOT EXISTS scp_idx ON scp_data (scp_name, scp_acc, scp_id)"
		self.curs.execute(sql)
		
		sql = "CREATE INDEX IF NOT EXISTS gak_index ON genome_accession_counts (genome_id, scp_id)"
		self.curs.execute(sql)
		
		self.conn.commit()
		
		if gen_type == "prokaryote":
			hmms, hmm_names = read_prok_hmm_directory(scp_dir)
		else:
			#scp file is the arg, but this is really a directory
			hmms, hmm_names = read_euk_hmm_directory(scp_dir)

			
		mn = hmm_jsonifier()
		mn.set_strings(hmms, hmm_names, gen_type)
		mn.jsonify()
		mn.format()
		json_fmt_hmms = mn.insertable_jsons
		
		sql = "INSERT OR IGNORE INTO scp_data VALUES (?, ?, ?, ?, ?)"
		self.curs.executemany(sql, json_fmt_hmms)
		self.conn.commit()
		
		
		
	def create_scp_table_pair(self, scp_label):
		genome_first = ', '.join([
				"genome_id INTEGER PRIMARY KEY",
				"tetramers BLOB"
				])
				
		tetra_first = ', '.join([
				"tetramer INTEGER PRIMARY KEY",
				"genomes BLOB"
				])
				
		sql = 'CREATE TABLE IF NOT EXISTS "{scp}_genomes" ({fields})'.format(scp = scp_label, fields = genome_first)
		self.curs.execute(sql)
		
		sql = 'CREATE TABLE IF NOT EXISTS "{scp}_tetras" ({fields})'.format(scp = scp_label, fields = tetra_first)
		self.curs.execute(sql)
		
		sql = 'CREATE INDEX IF NOT EXISTS "{scp}_genome_first" ON "{scp}_genomes" (genome_id)'.format(scp = scp_label)
		self.curs.execute(sql)
		
		sql = 'CREATE INDEX IF NOT EXISTS "{scp}_tetra_first" ON "{scp}_tetras" (tetramer)'.format(scp = scp_label)
		self.curs.execute(sql)
		
		
	def prep_scp_tables(self):
		sql = "SELECT scp_acc FROM scp_data"
		scp_names = self.curs.execute(sql).fetchall()
		for scp in scp_names:
			self.create_scp_table_pair(scp[0])
			print("\tSCP:", scp[0], "initialized!")
			
			
	def create_ruleset_table(self, rule_json = None, aai_models = None):
		#This is the set of rules to transform each SCP jacc index into AAI. 
		#Default transformation is multiply by 1 (e.g. unchanged Jacc. idx)
		rule_fields = ', '.join([
				"scp_name TEXT",
				"class1 INTEGER",
				"class2 INTEGER",
				"low_bound REAL",
				"high_bound REAL",
				"weight REAL",
				"intercept REAL",
				"jacc_slope REAL",
				"min_length_slope REAL",
				"length_diff_slope REAL",
				"min_score_slope REAL",
				"score_diff_slope REAL"
				])
		
		#Ruleset table for determining genome class in adding genomes to a database
		rule_addition_fields = ', '.join([
				'ruleset_json TEXT'
				])
				
		#Allows use of these values:
		#final_aai_or_jacc_estimate, genome_length1, genome_length2, hmm_score_1, hmm_score_2, shared_scp_ct, scp_std_dev
		#Currently all of these models must take the form of linear fits.
		post_processing_rules = ', '.join([
		'genome_class_1 INTEGER',
		'genome_class_2 INTEGER',
		'intercept REAL',
		'jacc_or_aai_slope REAL',
		'shared_scp_ct REAL',
		'scp_stddev REAL'])
		
		db_genomes_rule = None
		final_transform_rules = None
		if rule_json is None:
			default_rule = json.dumps({"class_assignment":None,
										"final_transform":None})
			rule = default_rule
		else:
			with open(rule_json) as fh:
				rule = fh.read()
			class_extraction = json.loads(rule)
			rule_update = []
			if 'class_assignment' in class_extraction:
				if 'classes' in class_extraction['class_assignment']:
					for c in class_extraction['class_assignment']['classes']:
						class_id = class_extraction['class_assignment']['classes'][c]["class"]
						rationale = class_extraction['class_assignment']['classes'][c]["rationale"]
						
						#Need to get the actual rule
						scp_membership = class_extraction['class_assignment']['classes'][c]['members']
						#transform to string
						scp_membership = '; '.join(scp_membership)
						
						next_rule = (class_id, c, scp_membership, rationale,)
						rule_update.append(next_rule)
						
				sql = "INSERT OR REPLACE INTO genome_classes VALUES (?, ?, ?, ?)"
				self.curs.executemany(sql, rule_update)
				self.conn.commit()
			
			if 'final_transform' in class_extraction:
				final_transform_rules = []
				for genome_class_1 in class_extraction['final_transform']:
					for genome_class_2 in class_extraction['final_transform'][genome_class_1]:
						this_dat = class_extraction['final_transform'][genome_class_1][genome_class_2]
						intercept = this_dat['intercept']
						jslope = this_dat['jacc_or_aai_slope']
						scp_ct_slope = this_dat['shared_scp_ct']
						scp_dev_slope = this_dat['scp_stddev']
						next_rule = (int(genome_class_1), int(genome_class_2), intercept, jslope,
									scp_ct_slope, scp_dev_slope,)
						final_transform_rules.append(next_rule)
				
		rule = (rule, )
		
		sql = "DROP TABLE IF EXISTS rule_addition_fields"
		self.curs.execute(sql)
		sql = "CREATE TABLE IF NOT EXISTS rule_addition_fields ({fields})".format(fields = rule_addition_fields)
		self.curs.execute(sql)
		sql = "INSERT INTO rule_addition_fields VALUES (?)"
		self.curs.execute(sql, rule)
		self.conn.commit()
		
		
		sql = "DROP TABLE IF EXISTS final_transform"
		self.curs.execute(sql)
		sql = "CREATE TABLE final_transform ({fields})".format(fields = post_processing_rules)
		self.curs.execute(sql)
		sql = "INSERT INTO final_transform VALUES (?, ?, ?, ?, ?, ?)"
		if final_transform_rules is not None:
			self.curs.executemany(sql, final_transform_rules)
			self.conn.commit()
		
		
		#Ruleset table can't have duplicates
		sql = "DROP TABLE IF EXISTS rulesets"
		self.curs.execute(sql)
		sql = "CREATE TABLE IF NOT EXISTS rulesets ({fields})".format(fields = rule_fields)
		self.curs.execute(sql)
		
		
		self.conn.commit()
		
		aai_transform_rules = []
		sql = "SELECT scp_acc FROM scp_data"
		recovered_scps = self.curs.execute(sql).fetchall()
		recovered_scps = [s[0] for s in recovered_scps]
		recovered_scps = set(recovered_scps)
		skip_default = False
		if aai_models is not None:
			skip_default = True
			detected_scps = []
			with open(aai_models) as fh:
				header = fh.readline()
				for line in fh:
					segs = line.strip().split("\t")
					scp = segs[0]
					detected_scps.append(scp)
					c1, c2 = int(segs[1]), int(segs[2]) #genome classes
					lb, hb = float(segs[3]), float(segs[4]) #low bound, high bound
					weight = float(segs[5])
					jacc_int, jacc_slope = float(segs[6]), float(segs[7])
					minlen, lendiff = float(segs[8]), float(segs[9])
					minscore, scorediff = float(segs[10]), float(segs[11])
					next_rule = (scp, c1, c2, lb, hb, weight, jacc_int, jacc_slope, minlen, lendiff, minscore, scorediff, )
					aai_transform_rules.append(next_rule)
			
			detected_scps = set(detected_scps)
			if detected_scps != recovered_scps:
				print("There was a difference in the set of supplied, per-SCP AAI transformation models and the HMMs used to create the database.")
				print("There must be a rule for EVERY SCP in the database if rules are to be manually added.")
				skip_default = False

		if not skip_default:
			#We cannot add an incomplete set of SCP rules
			aai_transform_rules = []
			for scp in recovered_scps:
				#Class zero is the default, low bound/high bound covers the entire range of possible jacc. values. 
				#slope and intercept are 1/0 so that the jacc. is functionally untransformed, other values zero to not affect results
				next_rule = (scp[0], #SCP label
							0, 0, #genome class 1, class 2
							0.0, 1.0, #low bound high bound on incoming jacc. range
							1.0, #model weight
							0.0, 1.0, #jaccard intercept, slope
							0.0, 0.0, #min length, length diff slopes
							0.0, 0.0,) #min score, score diff slopes
				aai_transform_rules.append(next_rule)
			
		sql = "INSERT INTO rulesets VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
		self.curs.executemany(sql, aai_transform_rules)
		self.conn.commit()
		
		sql = ''
		
		
	def prep_db(self, hmms, gtype = "prokaryote", rule_json = None, models = None):
		print("Initializing genome metadata")
		self.initialize_genome_index_table()
		print("Adding SCP records")
		self.initialize_gak(hmms, gtype)
		print("Adding default AAI calculation ruleset")
		self.create_ruleset_table(rule_json, models)
		print("Initializing SCP tables")
		self.prep_scp_tables()
		#Commit tables
		self.conn.commit()
		print("Done! FastAAI database created!")

import argparse		
def db_prep_options():
	#Needs
	#database_path -p / --path
	#hmm_files -m / --hmms
	#is_euk -e / --euk flag
	#rules_json -r / --rules
	#per_scp_aai_models -a / -aai_models
	parser = argparse.ArgumentParser(description='Options for initializing a FastAAI v2 database. Only --path and --hmms are required.')
	
	parser.add_argument('-d', '--path', dest = 'dbpath',
						help='A path to the FastAAI v2 database to be created by this command')
	parser.add_argument('-m', '--hmms', dest='hmms',
						help='Path to a directory containing the HMM models for each single copy protein this database will be initialized with.')
	parser.add_argument('-e', '--euk', dest='is_euk', action='store_true',
						help='Flag indicating whether the genomes are eukaryotic microbes.')
	parser.add_argument('-r', '--rules', dest='rules_json', default = None,
						help='A file containing rules for the assignment of genome classes to incoming genomes and post-processing on genomes to enhance AAI estimates')
	parser.add_argument('-a', '--aai_models', dest='aai_model_file', default = None,
						help='Tab-sep file containing piecewise linear models representing the per-SCP jaccard -> AAI transformation. May be replaced in the future with a more sophisticated model system.')

	args = parser.parse_known_args()
	return parser, vars(args[0])


def initialize_fastaai_database():
	parser, args = db_prep_options()
	if len(sys.argv) < 3:
		parser.print_help()
	else:
		f = args['dbpath']
		hmm_dir = args['hmms']
		is_euk = args['is_euk']
		rulej = args['rules_json']
		aai_models = args['aai_model_file']
		if is_euk:
			gtype = "eukaryote"
		else:
			gtype = "prokaryote"
		
		mn = fastaai_db_prepper(f)
		mn.open()
		mn.prep_db(hmm_dir, gtype, rulej, aai_models)
		mn.close()

