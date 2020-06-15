import sys
gff3 = open(sys.argv[1]+"/combined_gff3.tsv", "r+")
xgff3 = open(sys.argv[1]+"/combined_xgff3.tsv", "w+")

# ['Tzi8', 'Zm00042a019129', 'exon', '238300', '238927', '1\n']
# ['CML228', '632071', '633682\n']
for g in gff3:
    g = str(g).split("\t")
    xlim = open(sys.argv[1]+"/assigned_xlims.csv", "r+")
    for x in xlim:
        x = str(x).split(",")
        if g[0]==x[0] and g[2]=="mRNA":
            gl = int(g[3])
            gr = int(g[4])
            xl = int(x[1])
            xr = int(x[2][:-1]) # shouldn't it be -2?

            if gl<=xl and gr>=xr:
                xgff3.write("\t".join([g[0],g[1],g[2],str(xl),str(xr),g[5]]))
            if gl<=xl and gr>=xl and gr<=xr:
                xgff3.write("\t".join([g[0],g[1],g[2],str(xl),str(gr),g[5]]))
            if gl>=xl and gl<=xr and gr>=xr:
                xgff3.write("\t".join([g[0],g[1],g[2],str(gl),str(xr),g[5]]))
