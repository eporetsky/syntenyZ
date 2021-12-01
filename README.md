# syntenyZ
Custom scripts to make comparative genomic plots using genoPlotR

# Shiny app update
Just added version 0.2 of the syntenyZ shiny app.
It still requires some work but it is fully functional, even if error-prone.
There are 2 examples provided: Bx1 gene cluster and PAN1.
I will provide a direct link to the blastdb file once I upload them somewhere but FASTA file contains scripts to generate them all.
I will work on making a couple of slides for the Shiny app since it is not intuitive, and will try to make improvments to the app. 

# Pipeline steps
1. Blast query sequences against specified BLASTdbs
2. Filter the blast results to include only region of interest
3. Run the syntenyZ notebook to parse blast results and generate outputs
4. Run syntenyZ.R to generate the genoPlotR comparative figure
5. Edit PDF (with InkScape) to generate the figure for publication

# TODO
* Write a manual for using syntenyZ
* Create a template step-by-step example file
* Add detailed comments to code
* Many more...

# Contact
If you have any questions please feel free to contact me by mail or twitter @externelly.\
Contact me at eporetsky at ucsd.edu and I will try to help go through it.\
I will continue working on it when possible to make it more accesible.

# Examples
Taken from supplementary figure of:
\
Ding, Y., Weckwerth, P.R., Poretsky, E., Murphy, K.M., Sims, J., Saldivar, E., Christensen, S.A., Char, S.N., Yang, B., Tong, A., Shen, Z., Kremling, K.A., Buckler, E.S., Kono, T., Nelson, D.R., Bohlmann, J., Bakker, M.G., Vaughan, M.M., Khalil, A.S., Betsiashvili, M., Dressano, K., Köllner, T.G., Briggs, S.P., Zerbe, P., Schmelz, E.A. and Huffaker, A. (2020) Genetic elucidation of interconnected antibiotic pathways mediating maize innate immunity. Nat. Plants, 6, 1375–1388.
\
![Zx1-4](https://github.com/eporetsky/syntenyZ/blob/master/Tutorial/Zx1-4.png?raw=true)
\\

Taken from supplementary figure of:
\
Poretsky, E., Dressano, K., Weckwerth, P., Ruiz, M., Char, S.N., Shi, D., Abagyan, R., Yang, B. and Huffaker, A. (2020) Differential activities of maize plant elicitor peptides as mediators of immune signaling and herbivore resistance. Plant J, tpj.15022.
\
![ZmPROPEPs](https://github.com/eporetsky/syntenyZ/blob/master/Tutorial/ZmPROPEPs.png?raw=true)

