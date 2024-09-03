import json

class hmm_jsonifier:
	def __init__(self):
		self.models_dict = None
		self.acc_to_names = None
		self.models_jsons = None
		
		self.gentype = None
		
		self.insertable_jsons = None
		self.next_hmm_index = 0
		
	def set_strings(self, strdict, names, genome_type = "prokaryote"):
		self.models_dict = strdict
		self.acc_to_names = names
		self.gentype = genome_type
		
	def jsonify(self):
		for hmm in self.models_dict:
			self.models_dict[hmm] = json.dumps(self.models_dict[hmm], indent = 4)
			
	def format(self):
		self.insertable_jsons = []
		for hmm in self.models_dict:
			'''
			"scp_name TEXT",
			"scp_acc TEXT",
			"organism_type TEXT"
			"scp_id INTEGER",
			"scp_model TEXT"
			'''
			
			name = self.acc_to_names[hmm]
			
			next_json = (name, hmm, self.gentype, self.next_hmm_index, self.models_dict[hmm],)
			self.insertable_jsons.append(next_json)
			self.next_hmm_index += 1
				
	