# rmf 4.19.2018

# Purpose: remove primer and vector sequences from Sanger sequencing reads
# Goal: allow further processing on these files

# step 1: create reference fasta file containing primer and vector sequences to remove from reads
# primers_SangerSeq.fasta
# pCR4_vectorSeq.fasta

# step 2: write bbduk command to remove these sequences
#runBBDuk.sh

# step 3: run bbduk
SGE_Batch -c 'runBBDuk.sh' -r logBBDukTrim

# step 4: check efficacy of primer and vector sequence removal
