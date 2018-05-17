# rmf 5.17.2018

import sys, os
from Bio import SeqIO

# USAGE: 
usage = 'python ' + os.path.basename(__file__) + ' <fasta file> '
if len(sys.argv) != 2 or '-h' in sys.argv or '--help' in sys.argv:
    print usage
    sys.exit()

# SUBROUTINES

# this subroutine renames the defline IDs in the fasta file
def renameFastaDeflines(fastaFile):
    # extract each read's defline ID and sequence and store in SeqIO object 'fastaRead'
    fastaRead = SeqIO.parse(fastaFile,'fasta') # specify file name (fastaFile) and file type ('fasta')
    # access the defline ID stored in the 'fastaRead' object
    print fastaRead.id()

# ARGUMENTS and MAIN

# assign variable name to user input fasta file
fastaFile = sys.argv[1]

# call this subroutine, above
renameFastaDeflines(fastaFile)
