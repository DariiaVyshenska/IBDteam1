# rmf 5.17.2018
#!/local/cluster/bin/python
import sys, os  # for reading in files and running usage statement
import re  # for regular expressions, used to get taxonomic info out of BLAST files
from Bio import SeqIO # the following imports are for using various parts of the biopython package
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
    rec = []
    blastFile = open(blastFile)
    blastRecords = NCBIXML.parse(blastFile) # use parse() for multiple query sequences 
    for blastRecord in blastRecords:
        for alignment in blastRecord.alignments:
            line = alignment.title # assigns the alignment title containing taxonimic info to a variable
            rec.append(line) # append list of organisms
            break # only take the first one, break after that
    return rec

def extractNames(recs, regex):
    i = 0
    names = []
    for rec in recs:
        if regex.search(rec) != None:
            names.append(regex.search(rec).group(0))
        else:
            print("Genus not known", i)
        i += 1
    return names

# this subroutine renames the defline IDs in the fasta file
def renameFastaDeflines(organismNames,fastaFile):
    readList = []
    # set variable count to zero
    count = 0
    # store each read's defline ID and sequence in SeqIO object 'fastaReads'
    fastaReads = SeqIO.parse(fastaFile,'fasta') # specify file name (fastaFile) and file type ('fasta')
    # for each read in the object fastaReads
    for read in fastaReads:
        # GenBank defline format: Seq1 [organism=Genus species] [other info] description of the sequence
        newID = "Seq" + str(count) + " [organism=" + organismNames[count] + "]" # this variable must be a string
        # change the read ID
        read.name = ''
        read.description = ''
        read.id = newID
        #read = SeqRecord(read, id=newID) # newID variable must be a string in order to assign it as the record id
        readList.append(read)
        count +=1 # update count in preparation for next read
    return readList

def writeNewFasta(reads,outFile):
    with open(outFile, 'w') as outFile:
        SeqIO.write(reads,outFile,"fasta")

# ARGUMENTS and MAIN

# assign variable name to user input fasta file
fastaFile = sys.argv[1]
blastFile = sys.argv[2]

fastaFile_noPath = fastaFile.split('/')[-1]
fastaBase = fastaFile_noPath.split('.')[0]
outFile = str(fastaBase) + '_formatted.fasta'

regex = re.compile("[0-9a-zA-Z ]*;[0-9a-zA-Z -]*;[0-9a-zA-Z ._-]*$")

# call this subroutine, above
rec = parseBlastFile(blastFile)
names = extractNames(rec,regex) 
reads = renameFastaDeflines(names,fastaFile)
writeNewFasta(reads,outFile)
