#!/usr/bin/python3.7 -i
import torch
import esm
import sys
import subprocess
import os
import requests
import re
import time
import argparse
import biotite.structure.io as bsio
from Bio import SeqIO

my_parser = argparse.ArgumentParser(description='Benchmark ESMFold prediction with RMSD and TMScore metric for a known protein structure. Currently only supports pdb file format.')
my_parser.add_argument('-s', type=str, required=True, help = 'The known structure you want to benchmark ESM-Fold prediction capabilities')
my_parser.add_argument('-o', type=str, required=True, help = 'The output path where you want to save ESM-Fold predictions')
args = my_parser.parse_args()

known_structure = str(args.s)
out_path = str(args.o)
# autocompletion
import readline
import rlcompleter
readline.parse_and_bind('tab: complete')

# pymol launching
import pymol
pymol.pymol_argv = ['pymol','-qc']
pymol.finish_launching()
cmd = pymol.cmd

def esm_fold_for_large_seq(sequence):
	model = esm.pretrained.esmfold_v1()
	model = model.eval().cuda()
	with torch.no_grad():
		output = model.infer_pdb(sequence)
	return output

def get_tmscore_from_zang(native_structure,decoy_structure,path_to_tmscore_binary):
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
    cmd.load(str(native_structure), str(pdb_name))
    fasta_with_carrot = str(cmd.get_fastastr(str(pdb_name))).split("\n")
    fasta_without_carrot = "".join(fasta_with_carrot[1:])
    seq_length = len(fasta_without_carrot)
    predicted_structure_filename = path_to_save+pdb_name+'_esm_pred.pdb'
    pred_file_handler = open(predicted_structure_filename,"w+")
    pdb_string = esm_fold_for_large_seq(fasta_without_carrot)	
    pred_file_handler.write(pdb_string)
    pred_file_handler.close()
    is_pred_file_empty = os.stat(str(predicted_structure_filename)).st_size
    assert is_pred_file_empty != 0, "The prediction did not go through and the file is empty. Check script"
    decoy_existence = os.path.isfile(str(predicted_structure_filename))
    assert decoy_existence == True, "The prediction does not exist in the specified path."
    cmd.load(predicted_structure_filename,"esm_decoy")
    ca_rmsd  = round(float(str(cmd.align(str(pdb_name)+' and n. CA','esm_decoy and n. CA')).split()[3].split(',')[0]),3)
    tm_score = round(get_tmscore_from_zang(native_structure,predicted_structure_filename,"./TMscore"),3)
    plddt_column = bsio.load_structure(predicted_structure_filename, extra_fields=["b_factor"])
    mean_plddt = round(plddt_column.b_factor.mean(),3)
    print(native_structure,seq_length,str(ca_rmsd),str(tm_score),str(mean_plddt))
    cmd.quit()

benchmark_esmfold(known_structure,out_path)
