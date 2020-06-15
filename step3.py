import os, sys, pathlib
import pandas as pd


loci_list = os.listdir(sys.argv[1]+"/custom_loci")
loci_count = {}
for fl in loci_list:
	if fl==".DS_Store":
		continue
	fl_name = fl[:-4].split("_")[2]
	with open(sys.argv[1]+"/custom_loci/"+fl,"r") as loci_file:
		loci_count[fl_name] = 0
		for ln in loci_file:
			if ln.split(",")[2] == "mRNA":
				loci_count[fl_name] += 1

with open(sys.argv[1]+"/assigned_xlims.csv", "w") as assigned_xlims:
	assigned_xlims.write('genotype,left,right\n')
	for fl in loci_list:
		if fl==".DS_Store":
			continue
		fl_name = fl[:-4].split("_")[2]
		with open(sys.argv[1]+"/custom_loci/"+fl,"r") as loci_file:
			count = 0
			for ln in loci_file:
				if ln.split(",")[2] == "mRNA":
					assigned_xlims.write(fl_name+","+ln.split(",")[3]+","+ln.split(",")[4]+"\n")
					count += 1
