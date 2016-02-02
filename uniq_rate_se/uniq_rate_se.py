import sys
import getopt
import gzip
from re import sub
from Bio import SeqIO

def logging(*objs):
  print(*objs, file = sys.stderr)
  
  

def main(argv):

  try:
    opts, args = getopt.getopt(argv,"hi:",["ifile="])
  except getopt.GetoptError:
    logging('test.py -i <fasta.gz.file>')
    sys.exit(2)
  for opt, arg in opts:
    if opt == '-h':
      logging('test.py -i <fasta.gz.file>')
      sys.exit()
    elif opt in ("-i", "--ifile"):
      logging( 'Input file is "', arg)
      inputfile = gzip.open(arg, 'rt')
 
  outputfile = sys.stdout

  count = 0 
  uniq = 0
  seq_map = {}
  fasta_sequences = SeqIO.parse(inputfile ,'fastq')
  for fasta in fasta_sequences:
    count += 1
    if fasta.seq not in seq_map:
      uniq += 1
      seq_map[fasta.seq] = 0
  inputfile.close()
  print("Unique: ", uniq, "\nTotal: ", count)
  
      
   
if __name__ == "__main__":
   main(sys.argv[1:])
