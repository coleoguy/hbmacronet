Sample operations/commands of pipeline procedures 
procedure starts with a tree of 5 species with hybridizations already in the ms format via R
TODO write the first steps to get to a tree in ms format
overall file output order:
sample_ms -> gene_trees_sample -> dna_sample -> dna_sample.nex -> MB_output



## generating 4 gene trees by using a tree in ms format
./ms 5 4 -I 5 1 1 1 1 1 -ej 1.0 5 4 -ej 2.0 4 3 -ej 5.0 3 2 -em 5.0 1 2 250.0 -ej 10.0 1 2 -T | tail -n +4 | grep -v // > gene_trees_sample.txt

*go to output text file and add "[200]" on the front of each gene tree by hand -this is used in seq-gen* 
TODO write code code to this so it can be done for many times 



## Use gene trees to simulate 200bp per gene for each species
 ./seq-gen -mHKY -l 1000 -s 0.01 -p 5 < gene_trees_sample.txt > dna_sample.txt


*Convert output from seq-gen to Nexus type by hand and with aid of tool on http://sequenceconversion.bugaco.com/converter/biology/sequences/fasta_to_nexus.php
TODO Find a better way to automate this

##Commands for Mr Bayes to 
execute dna_sample.nex.txt
lset nst=6 rates=invgamma
mcmc ngen=20000 samplefreq=100 printfreq=100 diagnfreq=1000
sumt
