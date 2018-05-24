# rmf 5.17.2018

import sys, os  # for reading in files and running usage statement
import re  # for regular expressions, used to get taxonomic info out of BLAST files
# the following imports are for using various parts of the biopython package
from Bio import SeqIO  
from Bio.Blast import NCBIXML
from Bio.SeqRecord import SeqRecord

# USAGE: parses blast xml file and uses seq IDs to rename fasta file deflines
usage = 'python ' + os.path.basename(__file__) + ' <fasta file to modify> <BLAST XML outfile> \nUses sequence IDs from BLAST XML output file to modify FASTA file deflines in NCBI GenBank accepted format'
if len(sys.argv) != 3 or '-h' in sys.argv or '--help' in sys.argv:
    print usage
    sys.exit()

# SUBROUTINES

# parses blast xml output file
def parseBlastFile(blastFile):
    #organismName = re.compile('[;\s]([\D]+)[;\s]') # should get all taxonomic info
    blastFile = open(blastFile)
    blastRecords = NCBIXML.parse(blastFile) # use parse() for multiple query sequences 
    for blastRecord in blastRecords:
        for alignment in blastRecord.alignments:
            for hsp in alignment.hsps:
                line = alignment.title
                testList = line.strip().split(';')

# this subroutine renames the defline IDs in the fasta file
def renameFastaDeflines(fastaFile):
    # set variable count to zero
    count = 0
    # store each read's defline ID and sequence in SeqIO object 'fastaReads'
    fastaReads = SeqIO.parse(fastaFile,'fasta') # specify file name (fastaFile) and file type ('fasta')
    # for each read in the object fastaReads
    for read in fastaReads:
        # add one to count for each read in order to get unique number for each read
        count += 1
        # GenBank defline format: Seq1 [organism=Genus species] [other info] description of the sequence
        newID = "Seq" + str(count) + " [organism=" + organismName + "]" # this variable must be a string
        print read.id 
        # change the read ID
        read = SeqRecord(read, id=newID) # newID variable must be a string in order to assign it as the record id
        print read.id

# ARGUMENTS and MAIN

# assign variable name to user input fasta file
fastaFile = sys.argv[1]
blastFile = sys.argv[2]

# call this subroutine, above
parseBlastFile(blastFile)
renameFastaDeflines(fastaFile)
