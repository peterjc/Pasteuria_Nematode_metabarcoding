cd /home/jo42324/metabarcode/analysis/sequence-analysis/final-pipe/1_QC-and-assembly


# Generate a QC report using FastQC for paired end Illumina read data 
fastqc --nogroup /home/jo42324/metabarcode/data/raw_data/Raw_illumina/10623_metagen_f_reads.fastq -o trimqc/
wait
fastqc --nogroup home/jo42324/metabarcode/data/raw_data/Raw_illumina/10623_metagen_r_reads.fastq -o trimqc/


