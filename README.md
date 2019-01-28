# allelic_imbalance_analysis
# My project at the Furey Lab has been to study allelic imbalance in coding and noncoding regions in Crohn's Disease patient samples. A full description of the project is in the Allelic Imbalance Paper. 
# The Allelic Imbalance Pipeline image shows where each script is used
The create_personal_genome script is part of a process to convert a haploid reference genome to a custom diploid genome for each sample.
The allele_flipper script iterates through a .bam file, generating reads from the 5' to 3' direction that are outputted in one .fastq file, and other reads from 3' to 5' in another .fastq file. Each read has two copies; one with the reference allele and one with the alternate allele. These alleles and their locations are found in the other input, a .bed file. 
