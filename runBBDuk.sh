for f in /ACTF/Course/mb599_bds_s18/data/share/group1_Ream/IBDfasta/*; do infile=$(basename -- "$f"); file="${infile%.*}"; outfile="${file}.trimmed.fasta"; bbduk.sh in=$infile out=$outfile ref=pCR4_vectorSeq.fasta stats=sangerStats_vectorR.txt rcomp=t ktrim=r mink=11 k=31 hdist=1; done

for f in /ACTF/Course/mb599_bds_s18/data/share/group1_Ream/IBDfasta/*; do infile=$(basename -- "$f"); file="${infile%.*}"; outfile="${file}.trimmed.fasta"; bbduk.sh in=$infile out=$outfile ref=pCR4_vectorSeq.fasta stats=sangerStats_vectorL.txt rcomp=t ktrim=l mink=11 k=31 hdist=1; done

for f in /ACTF/Course/mb599_bds_s18/data/share/group1_Ream/IBDfasta/*; do infile=$(basename -- "$f"); file="${infile%.*}"; outfile="${file}.trimmed.fasta"; bbduk.sh in=$infile out=$outfile ref=primers_SangerSeq.fasta stats=sangerStats_primers.txt rcomp=t ktrim=r mink=11 k=20 hdist=1; done

# remove vector sequences on either side of insertion site (range 11 to 31 nt) from 3' end of reads
# bbduk.sh in=$infile out=$outfile ref=pCR4_vectorSeq.fasta stats=sangerStats_vectorR.txt rcomp=t ktrim=r mink=11 k=31 hdist=1

# remove vector sequences on either side of insertion site (range 11 to 31 nt) from 5' ends of reads 
# bbduk.sh in=$infile out=$outfile ref=pCR4_vectorSeq.fasta stats=sangerStats_vectorL.txt rcomp=t ktrim=l mink=11 k=31 hdist=1

# remove primers (range 11 to 20 nt) from 3' ends of reads
# bbduk.sh in=$infile out=$outfile ref=primers_SangerSeq.fasta stats=sangerStats_primers.txt rcomp=t ktrim=r mink=11 k=20 hdist=1

# rmf 4.19.2018
