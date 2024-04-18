#! /usr/bin/env python3

### This is a wrapper script to run several different tools on a set of 

import os

working_dir = '/home/piotr.lukasik/20240415_Wolbachia_TE'

path_to_raw_data = '%s/raw_data' % working_dir
path_to_filtered_data = '%s/filtered_data' % working_dir
path_to_reference_genomes = '%s/ref_genomes' % working_dir
   # have indexed them all - "for file in *fasta; do bwa index -a is $file; done"
   
path_to_qc_folder = '%s/qc_info' % working_dir

genome_list = [
['MG1', 'MG1-Dmauri_R1.fq.gz', 'MG1-Dmauri_R2.fq.gz', 'Wol_Dmauritiana.fasta'],
[['MG4', 'MG4-Dmelanogaster_R1.fq.gz', 'MG4-Dmelanogaster_R2.fq.gz', 'Wol_Dmelanogaster.fasta'],
['MG5', 'MG5-Brugiapahangi_R1.fq.gz', 'MG5-Brugiapahangi_R2.fq.gz', 'Wol_Brugiapahangi.fasta'],
['MG6', 'MG6-Armadilliumvulgare_R1.fq.gz', 'MG6-Armadilliumvulgare_R2.fq.gz', 'Wol_Avulgare.fasta'],
['TE1', 'TE1-Dmauri_R1.fq.gz', 'TE1-Dmauri_R2.fq.gz', 'Wol_Dmauritiana.fasta'],
['TE4', 'TE4-Dmelanogaster_R1.fq.gz', 'TE4-Dmelanogaster_R2.fq.gz', 'Wol_Dmelanogaster.fasta'],
['TE5', 'TE5-Brugiapahangi_R1.fq.gz', 'TE5-Brugiapahangi_R2.fq.gz', 'Wol_Brugiapahangi.fasta'],
['TE6', 'TE6-Armadilliumvulgare_R1.fq.gz', 'TE6-Armadilliumvulgare_R2.fq.gz', 'Wol_Avulgare.fasta']]

for genome in genome_list:
    print("Processing genome %s ........." % genome[0])

    ### pre-specifying sample names to facilitate work downstream ...
    ref_fasta = path_to_reference_genomes + "/" + genome[3]
    R1 = "%s/%s" % (path_to_raw_data, genome[1])
    R2 = "%s/%s" % (path_to_raw_data, genome[2])
    R1_trimmed_name = path_to_filtered_data + "/" + genome[1].split(".")[0] + "_val_1.fq.gz"
    R2_trimmed_name = path_to_filtered_data + "/"  + genome[2].split(".")[0] + "_val_2.fq.gz"
    R1_simpler_name = path_to_filtered_data + "/"  + genome[0] + "_R1.fq.gz"
    R2_simpler_name = path_to_filtered_data + "/"  + genome[0] + "_R2.fq.gz"   
    core_output_name =  path_to_qc_folder + "/" + genome[0]
    
    
    ### Trimming reads. Requiring high quality (cutoff 30), but accepting that they may be relatively short after trimming (100)
    print("Trimming reads...")
    print("trim_galore -q 30 -e 0.1 --cores 8 --length 100 -o %s --paired %s %s" % (path_to_filtered_data, R1, R2))
    os.system("trim_galore -q 30 -e 0.1 --cores 8 --length 100 -o %s --paired %s %s" % (path_to_filtered_data, R1, R2))
    os.system("mv %s %s && mv %s %s " % (R1_trimmed_name, R1_simpler_name, R2_trimmed_name, R2_simpler_name))
    
    
    ### Mapping reads to reference genomes
    print("Mapping reads...")
    os.system("bwa index -a is %s" % ref_fasta)
    os.system("bwa mem -t 40 -k 40 -T 80 %s %s %s | samtools view -b -S -F 12 - | samtools sort - > %s.bam" % (ref_fasta, R1_simpler_name, R2_simpler_name, core_output_name))
        # -k INT		minimum seed length [19] 
        # -A INT		score for a sequence match, which scales options -TdBOELU unless overridden [1]
        # -B INT		penalty for a mismatch [4]
        # -O INT[,INT]  gap open penalties for deletions and insertions [6,6]
        # -E INT[,INT]  gap extension penalty; a gap of size k cost '{-O} + {-E}*k' [1,1]
        # -L INT[,INT]  penalty for 5'- and 3'-end clipping [5,5]
        # -U INT        penalty for an unpaired read pair [17]
        # -T INT        minimum score to output [30]
        # -t INT        number of threads [1]
            
    ## Outputting mapping statistics
    print("Outputting stats...")    
    os.system("samtools index %s.bam" % core_output_name)
    os.system("samtools idxstats %s.bam > %s.info" % (core_output_name,core_output_name))
    
    ## Running qualimap
    print("Running Qualimap...")
    os.system("/home/piotr.lukasik/software/qualimap_v2.3/qualimap bamqc -c -ip -nt 60 -nw 1000 -bam %s.bam -outdir %s  -oc %s.cov -outformat pdf" % (core_output_name, path_to_qc_folder, core_output_name))
    os.system("mv %s/genome_results.txt %s/%s_results.txt" % (path_to_qc_folder, path_to_qc_folder, genome[0]))
    os.system("mv %s/report.pdf %s/%s_report.pdf" % (path_to_qc_folder, path_to_qc_folder, genome[0]))

    
    print("Genome %s successfully processed!" % genome[0])