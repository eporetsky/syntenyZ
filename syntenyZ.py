import os, sys

"""
sys.argv[1] - name of folder containing the blast hits
"""
if "1" in sys.argv[2]:
    print("turned step1 off for now")
#    os.system("python3 step1.py " + sys.argv[1])


"""
Step 2:
Input:
Output:
    1. coordinate_list.tsv
    2. custom_loci/ folder
    3. custom_genes/ folder
"""
if "2" in sys.argv[2]:
    os.system("python3 step2.py " + sys.argv[1])

"""
Step 3:
Input: custom_loci/ csv files
Output: assigned_xlims.csv
"""
if "3" in sys.argv[2]:
    os.system("python3 step3.py " + sys.argv[1])

"""
Step 4:
Input: coordinate_list.tsv
Output: a single fasta file for each locus
"""
if "4" in sys.argv[2]:
    print("turned step4 off for now")
    #os.system("python3 step4.py " + sys.argv[1])

"""
Step 5:
Input:
    1. coordinate_list.tsv
    2. assigned_xlims.csv
Output:
    1. combined_gff3.tsv
    2. combined_xgff3.tsv
Uses the edited gff3 tables in the GFF3 folder
"""
if "5" in sys.argv[2]:
    os.system("python3 step5.py " + sys.argv[1])

if "6" in sys.argv[2]:
    os.system("python3 step6.py " + sys.argv[1])
