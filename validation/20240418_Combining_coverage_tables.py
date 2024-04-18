#! /usr/bin/env python3

import os

   # The working directory on computing cluster
working_dir_on_azor = '/home/piotr.lukasik/20240415_Wolbachia_TE'
   # This is where reference genomes for read mapping are
path_to_reference_genomes = '%s/ref_genomes' % working_dir_on_azor
   ### Note: have indexed them all - "for file in *fasta; do bwa index -a is $file; done"
   
   # The primary data output folder
path_to_qc_folder = '%s/qc_info'  % working_dir_on_azor
   # Table with information on start and end position of each targetted region within each of the genomes
gene_pos_in_genomes = '%s/targets/gene_pos_in_ref_genomes.txt'  % working_dir_on_azor
   # The final output directory
path_to_vis = '/%s/qc_info/vis'  % working_dir_on_azor

genome_list = [
['Wol_Dmauritiana', 'MG1', 'TE1'],
['Wol_Dmelanogaster', 'MG4', 'TE4'],
['Wol_Brugiapahangi', 'MG5', 'TE5'],
['Wol_Avulgare', 'MG6', 'TE6']]

def ImportFasta(fasta_file):
   FASTA = open(fasta_file, 'r')
   Seq_list = []
   Sequence = ''
   Seq_heading = ''
   for line in FASTA:   # Copying the sequence (potentially spread across multiple lines) to a single line
      if line.startswith('>'):
         if Sequence != '':   # Saves the existing Seq_heading and Sequence to a list before overwriting them
            Seq_list.append([Seq_heading, Sequence])
         Sequence = ''
         Seq_heading = line.strip().strip(">").split()[0] # Takes everything until the first space as heading - stripping redundant info afterwards
         # This makes sense - here's an example heading: ">NODE_2_length_1018_cov_3.595960/910-218 [subseq from] NODE_2_length_1018_cov_3.595960"
      elif len(line.strip()) > 0:
         Sequence = Sequence + line.strip().upper()
   Seq_list.append([Seq_heading, Sequence]) # Saves the final sequence (Seq_heading and Sequence) to a list
   FASTA.close()
   return(Seq_list)
   
def ImportTable(table_file):
    TABLE = open(table_file, 'r')
    TAB = []
    for line in TABLE:
        if len(line.strip()) > 0:
            TAB.append(line.strip().split())
    return(TAB)
    





for genome in genome_list:
    
    print("Processing genome %s ..." % genome[0])
    
    #### Creating table (list-of-lists) where every base in any of the contigs is represented by a separate list
    #### Some genomes are represented as single contigs - easy - but others by multiple
    #### Hence, to keep track of the base_no corresponding to the start of each contig, I am creating a dictionary

    coverage_by_position = [['Contig', 'Pos']]
    contig_start_positions = {}
    
    # Initially, creating a table with contig names and base_nos
    print("Creating a table with contig names and base_nos...")
    
    FASTA = ImportFasta("%s/%s.fasta" % (path_to_reference_genomes, genome[0]))
    base_no = 0
    for seq in FASTA:
        
        ### Add info about the number of the first base in each contig to the "contig_start_positions" dictionary
        contig_start_positions[seq[0]] = base_no+1
        
        ### For every base in the assembly, add contig_name and base_no to the list
        for pos in range(len(seq[1])):
            base_no += 1
            coverage_by_position.append([seq[0], base_no])
            
    # printing summary of contigs/sized... for troubleshooting mostly!
    print("Table created; contig no = %s; total length = %s" % (len(contig_start_positions), len(coverage_by_position)-1))
    print("%s contigs found:" % len(contig_start_positions))
    
    for contig in contig_start_positions.keys():
        print("   %s starting at position %s" % (contig, contig_start_positions[contig]))
    print("")
        
    
    
    
    
    #### Adding info on bases spanned by target regions
    print("Adding info on bases spanned by target regions...")
    
    # Initially, just add xeros everywhere
    coverage_by_position[0].append("Target")
    for position in range(1,len(coverage_by_position)):
        coverage_by_position[position].append(0)
    
    # Reading table with info on coverage for any non-zero-coverage bases
    TARGET_positions = ImportTable(gene_pos_in_genomes)
    
    for line in TARGET_positions:
    ### lines look like this ---
    ### path_RefSeq_WP_006014985	Wol_Brugiapahangi	CP050521.1_FR3	267403	267630	227
        
        if line[1] == genome[0]:
            Wolb_contig = str(line[2].strip())
            if Wolb_contig != "0":         # 0 indicates regions not found in a given genome
                first_gene_coord = int(line[3].strip())
                second_gene_coord = int(line[4].strip())       
                gene_length = int(line[5].strip())
                                 
                gene_start_in_contig = min(first_gene_coord, second_gene_coord)
                gene_end_in_contig = max(first_gene_coord, second_gene_coord)
                
                target_name = str(line[0].strip())
                Wolb_genome = str(line[1].strip())
                
                
                position_start_corrected = int(contig_start_positions[line[2]])+gene_start_in_contig-1
                
                i = 0
                while i < (gene_length+1):
                    coverage_by_position[position_start_corrected+i][2] = line[0]
                    i += 1
            
            
            
    ### Now, positions are classified as 0 or "target", but I want to add level 1 - "neighbourhood", anything within 500bp from "target"

    print("Adding info on whether a position can be considered 'neighborhood'...")

    for pos in range(1, len(coverage_by_position)):
        if coverage_by_position[pos][2] == 0:
            start = max(1, pos-501)
            end = min(len(coverage_by_position), pos+501)
            for i in range(start, end, 50):
                if coverage_by_position[i][2] not in [0, 1] and coverage_by_position[pos][0] == coverage_by_position[i][0]:
                    coverage_by_position[pos][2] = 1
                    break
    




    ### Now, adding info on coverage in genomes:
    print("Adding info on coverage in genomes...")

    
    for lib_no in range(1,len(genome)):
        lib = genome[lib_no]  ### e.g., 'MG1'
        print("Adding info - genome %s..." % lib)
        
        coverage_by_position[0].append(lib)
        for position in range(1,len(coverage_by_position)):
            coverage_by_position[position].append(0)
        
        genome_coverage_table = ImportTable("%s/%s.cov" % (path_to_qc_folder, lib))
        for line in genome_coverage_table:
            # Except for the first line, goes like this: ""#chr pos coverage" -
            # NC_002978.6_wMel	8	172 
            
            if not line[0].startswith("#"):        
                contig_name = str(line[0])
                pos = int(line[1])+contig_start_positions[contig_name]-1
                cov = int(line[2])
                
                coverage_by_position[pos][-1] = cov
     
     
     
    
    #### Finally, exporting the table ... 
    print("Exporting full data table!")
    export_table = open("%s/%s_coverage_info.txt" % (path_to_vis, genome[0]), "w")

    for row_no in range(len(coverage_by_position)):
        for col_no in range(len(coverage_by_position[0])-1):
            print(coverage_by_position[row_no][col_no], end = "\t", file=export_table)
        print(coverage_by_position[row_no][-1], file=export_table)
    
    export_table.close()
    
    
    
    
    ### Now, compacting the table ... computing summary stats for every window of 50bp
    print("Compacting data table ...")
        
    vis_table = [coverage_by_position[0]]
    
    for i in range(1, len(coverage_by_position)-51, 50):
        line_from_coverage_by_position = coverage_by_position[i]
        # e.g., ['CP034335_wMau', 750000, 'phy_RefSeq_WP_006279887', 51, 4594]
        # or    ['CP034335_wMau', 950000, 1, 48, 0]
        contig_name = line_from_coverage_by_position[0]
        if coverage_by_position[i+50][0] == contig_name:
            #classify bin of 50: >1 target -> target; >1 "1" -> neighborhood; otherwise 0
            new_line = [contig_name, i]
            
            bin_classification = 0
            
            for x in range(0,51,10):
                if coverage_by_position[i+x][2] == 1:
                    bin_classification = 1
                    break
                elif coverage_by_position[i+x][2] not in [0, 1]:
                    bin_classification = coverage_by_position[i+x][2]
                    break
            
            new_line.append(bin_classification)
            
            for col in range(3,len(coverage_by_position[0])):
                cov_total_over_50bp = 0
                for x in range(0,50):
                    cov_total_over_50bp += coverage_by_position[i+x][col]
                cov_average = float(cov_total_over_50bp)/50
                    
                new_line.append(cov_average)
            
            vis_table.append(new_line)


    ### Exporting "compacted" data table!
    print("Exporting compacted data table!")

    export_compact_table = open("%s/%s_compacted_table.txt" % (path_to_vis, genome[0]), "w")

    for row_no in range(len(vis_table)):
        for col_no in range(len(vis_table[0])-1):
            print(vis_table[row_no][col_no], end = "\t", file=export_compact_table)
        print(vis_table[row_no][-1], file=export_compact_table)
    
    export_compact_table.close()
            
    
print("It's all done!")
