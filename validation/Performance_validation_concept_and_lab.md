# Performance validation

To validate the approach, we have conducted Wolbachia target enrichment on six "reference" samples, insects infected by different Wolbachia strains representing four supergroups that should
have had their genomes completely sequenced.
The strains are listed in [Table S1](Validation_Table_S1.txt).




### Navigate to the directory with alignments
```
cd /Users/piotrlukasik/bioinfo/TE_Wolbachia/targets
cd pathway_genes
```

### Create hmm profiles based on nuclotide alignments of marker genes
```
mkdir hmms
for file in *fa; do echo $file; hmmbuild --dna hmms/"$file".hmm $file; done
```

### Use hmm profiles to search your genomes and output the top hit (in my case, the genomes are in ~/bioinfo/TE_Wolbachia/ref_genomes/)
```
cd hmms
for file in *hmm; do for genome in ~/bioinfo/TE_Wolbachia/ref_genomes/*fasta; do echo "$file"___"$genome" | sed 's/\/Users\/piotrlukasik\/bioinfo\/TE_Wolbachia\/ref_genomes\///g' >> ../../pathways.hmmhits; nhmmer $file $genome | grep -A 3 "Scores for complete hits" | tail -1 >> ../../pathways.hmmhits; done; done
```

### Open files in your favorite text editor and use REGEX to change the format of files for comfortable viewing in Excel :)
```
find: .fasta\s+([ \t]+)(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+).*    replace: \t\1\t\2\t\3\t\4\t\5\t\6
find: .fasta\s+      replace: \n
find: .1.fa.hmm___   replace: \t
find: .fa.hmm___     replace: \t
```

### Repeat the process for phylogenetically informative and for special genes


