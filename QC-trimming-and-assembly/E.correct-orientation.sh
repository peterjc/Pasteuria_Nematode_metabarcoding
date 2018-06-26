cd /home/jo42324/metabarcode/analysis/sequence-analysis/Pasteuria_Nematode_metabarcoding/QC-trimming-and-assembly

# invoke the correct orientation python script
python correct-orientation.py
# wait for it to finish
wait 

# concatenate flipped reads and reads which assembled in the expected orientation 
cat pas-amplicons-orientation-corrected.fastq Pas-amplicons-containing-FandR-correct-orientation.fastq > all-pas-amplicons.fastq
wait 

# For both primer sets 
cat nem-amplicons-orientation-corrected.fastq Nem-amplicons-containing-FandR-correct-orientation.fastq > all-nem-amplicons.fastq
