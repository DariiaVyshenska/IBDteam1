for f in /ACTF/Course/mb599_bds_s18/data/share/group1_Ream/IBDfasta/*.fasta; do file=$(basename -- "$f"); filebase="${file%.*}"; outfile="bbdukTrim_out/${filebase}.trimmedR.fasta"; statsFile="bbdukTrim_stats/${filebase}_sangerStats_vectorR.txt"; bbduk.sh in=$f out=$outfile ref=pCR4_vectorSeq.fasta stats=$statsFile rcomp=t ktrim=r mink=11 k=31 hdist=1; done

for f in /ACTF/Course/mb599_bds_s18/home/feyr/IBDteam1/bbdukTrim/bbdukTrim_out/*.trimmedR.fasta; do file=$(basename -- "$f"); filebase="${file%.trimmedR.fasta}"; outfile="bbdukTrim_out/${filebase}.trimmedRL.fasta"; statsFile="bbdukTrim_stats/${filebase}_sangerStats_vectorL.txt"; bbduk.sh in=$f out=$outfile ref=pCR4_vectorSeq.fasta stats=$statsFile rcomp=t ktrim=l mink=11 k=31 hdist=1; done

for f in /ACTF/Course/mb599_bds_s18/home/feyr/IBDteam1/bbdukTrim/bbdukTrim_out/*.trimmedRL.fasta; do file=$(basename -- "$f"); filebase="${file%.trimmedRL.fasta}"; outfile="bbdukTrim_out/${filebase}.trimmedRLP.fasta"; statsFile="bbdukTrim_stats/${filebase}_sangerStats_primers.txt"; bbduk.sh in=$f out=$outfile ref=primers_SangerSeq.fasta stats=$statsFile rcomp=t ktrim=r mink=11 k=20 hdist=1; done

# remove vector sequences on either side of insertion site (range 11 to 31 nt) from 3' end of reads
# bbduk.sh in=$infile out=$outfile ref=pCR4_vectorSeq.fasta stats=sangerStats_vectorR.txt rcomp=t ktrim=r mink=11 k=31 hdist=1

# remove vector sequences on either side of insertion site (range 11 to 31 nt) from 5' ends of reads 
# bbduk.sh in=$infile out=$outfile ref=pCR4_vectorSeq.fasta stats=sangerStats_vectorL.txt rcomp=t ktrim=l mink=11 k=31 hdist=1

# remove primers (range 11 to 20 nt) from 3' ends of reads
# bbduk.sh in=$infile out=$outfile ref=primers_SangerSeq.fasta stats=sangerStats_primers.txt rcomp=t ktrim=r mink=11 k=20 hdist=1

# rmf 4.19.2018
