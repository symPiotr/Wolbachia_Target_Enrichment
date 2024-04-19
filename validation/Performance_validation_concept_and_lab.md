# Performance validation

### 1. Validation of the target enrichment process

To validate the approach, we have conducted Wolbachia target enrichment using DNA for four "reference" samples, invertebrates infected by different Wolbachia strains representing four supergroups that have had their genomes completely sequenced. 
The species are listed in [Table S1](Validation_Table_S1.txt).
For each of these samples, we prepared a library using NEBNext UltraExpress DNA kit, using half of the recommended reaction volume. After ligation of adapter stub and before indexing, the reaction was further divided into two halves, each of which was indexed with a different set of dual-unique indexes. Subsequently, for each sample, we used one of the resulting libraries for target enrichment reaction. For that, we prepared two pools, each including two experimental libraries and two others, combined at equimolar volumes. 
For the enrichment reaction, we followed the Standard protocol from MyBaits Manual v.5.0.3, with the hybridization temparature set to 60C.  
  
We combined the post-enrichment libraries with uncaptured libraries for the same samples for sequencing in an Illumina MiSeq v3 600-cycle flow cell. 

The reads were trimmed with [Trim_Galore v. 0.6.10](https://github.com/FelixKrueger/TrimGalore), with the quality cutoff of 30 and length cutoff of 100bp.
Then, we used captured and uncaptured reads from each library for mapping against their respective reference *Wolbachia* genomes, using [bwa v. 0.7.17-r1188](https://github.com/lh3/bwa), with seed length set to 40 and minimum outputted score of 80, which largely resolved the issue of non-specific read mapping to lower-complexity regions of the genome. We then used [samtools v. 1.18](https://www.htslib.org/) to filter unmapped reads and sort bam files, and [qualimap v. 2.3](http://qualimap.conesalab.org/) to assess per-base coverage of the genomes. 
The wrapper script for the steps above is provided - [20240418_mapping_reads_to_Wolbachia_genomes.py](20240418_mapping_reads_to_Wolbachia_genomes.py).  
  
We then combined information on *Wolbachia* genome coverage based on Target Enrichment and on unenriched Metagenomic dataset for the same biological samples with information on location of the targetted regions without the reference Wolbachia genomes. We then computed the average genomic coverage in 50-bp windows. Simultaneoulsy, we extracted information for each 50bp window whether any of its bases corresponds to the targetted region, or is within 500bp from the targetted region; if so, we classified the window as "on target", or "within target neighbourhood". We then visualized the coverage across the windows, within the initial 100kb of each genome, and for different target categories (phylogenetic, pathway, special) within entire genomes, using a custom Processing VERSION script. The Python script for combining information from across tables, the final tables for the four genomes, and the Processing script for data visualization are provided.  
  * [20240418_Combining_coverage_tables.py](ref)
  * [20240418_Genome_coverage_visualization.pyde](ref)
  * Data_tables: [Drosophila melanogaster / wMel](Wol_Dmelanogaster_compacted_table.txt), [Drosophila mauritiana / wXXX](Wol_Dmauritiana_compacted_table.txt), [Brugia pahangi / wXXX](Wol_Brugiapahangi_compacted_table.txt), [Armadillium vulgare / wXXX](Wol_Avulgare_compacted_table.txt).

The Combining_coverage_tables script also computes the average coverage for e
    
### 2. Testing the ability to reconstruct genes and pathways
  
For all targets, we have compared their presence within the reference genomes, and within the metagenomic and target enrichment datasets.
We used [nhmmer - part of HMMER 3.4](http://hmmer.org/) to search for each target within each genome, with e-value set to e-10, and used esl-sfetch to extract the target region coordinates and sequences.

We have used bam files constructed by mapping reads to the reference genomes, as explained above, then used samtools view to export reads mapped to each targetted region, and samtools fastq to export these reads.







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


