#Vsearch dereplication 

# Now that we have sorted and filtered our assembled reads (amplicons) we need to take all the amplicons binned (sorted into files for each primer pair) in the previous steps and de-replicate them at the study level (collapse all identical amplicons into a single record annotated with the total number of times it occurs). #

# change into the directory where your binned and filtered fastq files are 
cd /home/jo42324/metabarcode/analysis/sequence-analysis/final-pipe/1_QC-and-assembly
  
# invoke vsearch to dereplicate reads at the study level  
vsearch --derep_fulllength all-pas-amp-maxeefiltered.fasta \
--output global-pas-derep.fasta \
--sizeout 

wait

vsearch --derep_fulllength all-nem-amp-maxeefiltered.fasta \
--output global-nem-derep.fasta \
--sizeout 


# VSEARCH command line options #
#Dereplication and rereplication
#  --derep_fulllength FILENAME dereplicate sequences in the given FASTA file
#  --derep_prefix FILENAME     dereplicate sequences in file based on prefixes
#  --rereplicate FILENAME      rereplicate sequences in the given FASTA file
#Options
#  --maxuniquesize INT         maximum abundance for output from dereplication
#  --minuniquesize INT         minimum abundance for output from dereplication
#  --output FILENAME           output FASTA file
#  --relabel STRING            relabel with this prefix string
#  --relabel_keep              keep the old label after the new when relabelling
#  --relabel_md5               relabel with md5 digest of normalized sequence
#  --relabel_sha1              relabel with sha1 digest of normalized sequence
#  --sizein                    propagate abundance annotation from input
#  --sizeout                   write abundance annotation to output
#  --strand plus|both          dereplicate plus or both strands (plus)
#  --topn INT                  output only n most abundant sequences after derep
#  --uc FILENAME               filename for UCLUST-like dereplication output
#  --xsize                     strip abundance information in derep output