#!/usr/bin/python

import sys
import pathlib

project_name = sys.argv[1]
# print(project_name)

pathlib.Path("projects/"+project_name+"/").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_loci").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/fasta_loci").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_loci").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_CDS").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_CDS/genes").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_CDS/inbreds").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_intergenic").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_intergenic/genes").mkdir(parents=True, exist_ok=True)
pathlib.Path("projects/"+project_name+"/custom_intergenic/inbreds").mkdir(parents=True, exist_ok=True)
