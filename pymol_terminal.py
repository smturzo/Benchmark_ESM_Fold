#!/usr/bin/python3.7 -i

import sys
import subprocess
import os
import requests
import re
import time
import argparse
from Bio import SeqIO

#ABS_PATH = os.path.abspath('./')

# autocompletion
import readline
import rlcompleter
readline.parse_and_bind('tab: complete')

# pymol launching
import pymol
pymol.pymol_argv = ['pymol','-qc'] + sys.argv[1:]
pymol.finish_launching()
cmd = pymol.cmd
print("Printing Okay because we have successfully initialized pymol")
#cmd.load("./6DRV_A.pdb","test_file")
#my_fasta = str(cmd.get_fastastr("test_file")).split("\n")#.replace("\n","")
#my_fasta_string_only = "".join(my_fasta[1:])
#print(len(my_fasta_string_only))

def get_tmscore_from_zang(native_structure,decoy_structure,path_to_tmscore_binary="./TMscore"):
	assert path_to_tmscore_binary == True, "TMscore binary file not found in path. TMscore was not calculated and the script will now crash."
	tm_output = subprocess.check_output([str(path_to_tmscore_binary),str(native_structure),str(decoy_structure)])
	tm_output = str(tm_output).split("\\n")
	tm_score = None
	for d in tm_output:
		x = re.sub(r"\s\s+", " ", d).split(' ')
		if x[0] == "TM-score" and x[1] == "=":
			tm_score = float(x[2])

	return tm_score


def benchmark_esmfold(native_structure,path_to_save):
	native_existence = os.path.isfile(str(native_structure))
	assert native_existence == True, "Native structure does not exist in the specified path. Provide native structure with the full directory path"
	cmd.delete('all')
	pdb_name = str(native_structure).split("/")[-1].split(".pdb")[0]
	print(pdb_name)
	cmd.load(str(native_structure), str(pdb_name))
	fasta_with_carrot = str(cmd.get_fastastr(str(pdb_name))).split("\n")
	fasta_without_carrot = "".join(fasta_with_carrot[1:])
	print(len(fasta_without_carrot))
	headers = {
		'Content-Type': 'application/x-www-form-urlencoded',
	}

	response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=fasta_without_carrot)
	predicted_structure_filename = os.path.join(path_to_save, pdb_name) + "_esm_pred.pdb"
	pdb_string = response.content.decode('utf-8')
	if pdb_string.startswith("HEADER"):
		with open(predicted_structure_filename, "w") as out:
			out.write(pdb_string)
		print(f"Results saved to {predicted_structure_filename}")
	else:
		print(pdb_string)	
	
	decoy_existence = os.path.isfile(str(predicted_structure_filename))
	assert decoy_existence == True, "The prediction does not exist in the specified path. Dive into the script and let me know what I did wrong"
	cmd.load(predicted_structure_filename,"esm_decoy")
	ca_rmsd  = float(str(cmd.align(str(pdb_name)+' and n. CA','esm_decoy and n. CA')).split()[3].split(',')[0])
	tm_score = get_tmscore_from_zang(native_structure,predicted_structure_filename)
	print(native_structure,str(ca_rmsd),str(tm_score))

benchmark_esmfold("./2mlt_A.pdb","./")


