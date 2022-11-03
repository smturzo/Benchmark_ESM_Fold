#!/usr/bin/env python
from pymol import cmd
import os
import requests
import re
import time

ABS_PATH = os.path.abspath('./')

def query_esmfold(sequence:str, name:str=None):
    """Predict protein structure with ESMFold

    Args:
        sequence (str): amino acid sequence
        name (str, optional): _description_. Defaults to None.
    """
    sequence = re.sub("[^A-Z:]", "", sequence.replace("/",":").upper())
    sequence = re.sub(":+",":",sequence)
    sequence = re.sub("^[:]+","",sequence)
    sequence = re.sub("[:]+$","",sequence)
    print(sequence)
    headers = {
        'Content-Type': 'application/x-www-form-urlencoded',
    }


    response = requests.post('https://api.esmatlas.com/foldSequence/v1/pdb/', headers=headers, data=sequence)
    if not name:
        name = sequence[:3] + sequence[-3:] 
    pdb_filename = os.path.join(ABS_PATH, name) + ".pdb"
    pdb_string = response.content.decode('utf-8')
    if pdb_string.startswith("HEADER"):
        with open(pdb_filename, "w") as out:
            out.write(pdb_string)
        print(f"Results saved to {pdb_filename}")
    else:
        print(pdb_string)
    cmd.load(pdb_filename)
cmd.extend("esmfold", query_esmfold)
