from .agnostic_reader import agnostic_reader
import os

def read_prok_hmm_directory(hmmdir):
	hmm_contents = os.listdir(hmmdir)
	contents = []
	ar = None
	for f in hmm_contents:
		pp = os.path.normpath(hmmdir + "/" + f)
		ar = agnostic_reader(pp)
		next_file = ar.read()
		contents.append(next_file)
		
	contents = ''.join(contents)
	ar = None
	
	contents = contents.splitlines()
	
	current_model_acc = None
	scpname = None
	current_model_str = []
	acc_to_name = {}
	
	hmm_model_by_acc = {}
	for line in contents:
		if line.startswith("HMMER3"):
			if current_model_acc is not None:
				#print("A")
				current_model_str = '\n'.join(current_model_str)
				hmm_model_by_acc[current_model_acc] = current_model_str
				acc_to_name[current_model_acc] = scpname
			
			current_model_acc = None
			current_model_str = []
			
		if line.startswith("ACC"):
			accname = line.strip().split()[1:]
			accname = " ".join(accname)
			current_model_acc = accname
			
		if line.startswith("NAME"):
			scpname = line.strip().split()[1:]
			scpname = " ".join(scpname)
		
		current_model_str.append(line)
			
	#Final iteration	
	if current_model_acc is not None:
		#print("f")
		current_model_str = '\n'.join(current_model_str)
		hmm_model_by_acc[current_model_acc] = current_model_str
		acc_to_name[current_model_acc] = scpname
		
	return hmm_model_by_acc, acc_to_name

'''
arch_scps = ["TIGR03671",  "PF01994.16", "TIGR02338",  "PF01982.16", "PF04010.13", "TIGR01952",  "TIGR03626", 
			"TIGR02390",  "TIGR01052",  "PF01984.20", "TIGR03680",  "TIGR01028",  "TIGR03670",  "PF04104.14"
			"TIGR03674",  "PF00833.18", "TIGR00522",  "TIGR00491",  "TIGR01012",  "PF01667.17", "TIGR00483", 
			"PF01198.19", "PF01981.16", "PF01866.17", "TIGR01008",  "PF01194.17", "TIGR02236",  "PF00935.19"
			"TIGR01046",  "PF01282.19", "TIGR03676",  "PF06093.13", "TIGR01020",  "PF00900.21", "PF16906.5", 
			"TIGR00982",  "TIGR00037",  "PF02996.17", "TIGR00448",  "PF07541.13", "TIGR00134",  "PF01201.22"
			"PF17144.4",  "PF01015.19", "PF00827.18", "PF01655.18", "PF01280.21", "PF01780.19", "TIGR00291", 
			"PF01090.20", "TIGR00405",  "PF01200.19", "PF01157.18", "TIGR01018",  "PF04919.13"]
			
bact_scps = ["TIGR00095",  "TIGR00634",  "TIGR00416",  "TIGR00186",  "TIGR00138",  "PF02576.18", "TIGR02273", 
			"TIGR01128",  "TIGR00083",  "PF02699.15", "TIGR01039",  "TIGR00539",  "TIGR00436",  "TIGR01394", 
			"TIGR02386",  "TIGR00084",  "PF00213.18", "PF00825.18", "TIGR00580",  "PF00119.20", "TIGR00382", 
			"TIGR00959",  "TIGR01082",  "TIGR01169",  "TIGR03263",  "TIGR00054",  "TIGR03725",  "TIGR01029", 
			"TIGR00922",  "TIGR01044",  "TIGR02075",  "TIGR03632",  "TIGR02013",  "TIGR01021",  "TIGR00090", 
			"TIGR03953",  "TIGR01632",  "TIGR00020",  "TIGR01510",  "TIGR01079",  "TIGR01164",  "TIGR03625", 
			"TIGR01009",  "PF02565.15", "PF01783.23", "TIGR00115",  "TIGR03654",  "TIGR00092",  "TIGR00460", 
			"TIGR01071",  "PF03726.15", "TIGR02191",  "TIGR00445",  "TIGR00593",  "TIGR00362",  "TIGR00019", 
			"TIGR01951",  "TIGR03594",  "TIGR02432",  "TIGR00431",  "TIGR00963",  "PF03652.15", "TIGR01017", 
			"TIGR01087",  "TIGR00472",  "TIGR00615",  "TIGR00420",  "PF00830.19", "TIGR02027",  "TIGR00635", 
			"PF01632.19", "TIGR02729",  "TIGR01393",  "TIGR00059",  "TIGR00663",  "TIGR03723",  "TIGR00158", 
			"TIGR01066",  "TIGR01391",  "TIGR02397",  "TIGR02012",  "PF02130.17", "PF01668.18", "TIGR00168", 
			"TIGR00088",  "TIGR00116",  "TIGR01953",  "TIGR00487",  "TIGR00496",  "PF01649.18", "PF01250.17"
			"PF01016.19", "PF00829.21", "PF00886.19", "PF01245.20", "PF01715.17", "PF01195.19", "TIGR00006"]
			
			
cross_scps = ["PF00203.21", "PF00177.21", "PF00297.22", "PF00276.20", "PF00366.20", "PF00252.18", "PF00164.25"
			"PF00573.22", "PF00238.19", "PF00466.21", "PF00410.20", "PF00347.23", "PF01725.16", "PF00411.19"
			"PF00416.22", "PF00380.19", "PF00572.18", "PF00312.22", "PF00318.20", "PF00750.19", "TIGR00337", 
			"PF00338.22", "PF00281.19", "TIGR00065",  "TIGR00392",  "TIGR00468",  "PF01176.19", "TIGR00442", 
			"TIGR00398",  "TIGR00435",  "TIGR00456",  "TIGR00344",  "PF00162.19", "PF13393.6",  "PF01351.18"
			"PF01000.27", "PF00749.21"]



def filter(file, group, output):
	
	mod, an = read_hmm_model(file)

	with open(output, "w") as out:
		for scp in mod:
			if scp in group:
				print(mod[scp], file = out)
	

a = "models/finalized_archaeal_hmms.txt"
b = "models/finalized_bacterial_hmms.txt"
c = "models/finalized_cross_domain_hmms.txt"

filter(a, arch_scps, "filtered_to_candidates_hmm_models/archaeal_scps.hmm.txt")
filter(b, bact_scps, "filtered_to_candidates_hmm_models/bacterial_scps.hmm.txt")
filter(c, cross_scps, "filtered_to_candidates_hmm_models/cross_domain_scps.hmm.txt")
'''