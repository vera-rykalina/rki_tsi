#!/usr/bin/env python
from __future__ import print_function
from Bio import SeqIO, Seq
from re import sub
from AuxiliaryFunctions import ungap
import argparse
import os
import sys

## Author: Chris Wymant, chris.wymant@bdi.ox.ac.uk
## Acknowledgement: I wrote this while funded by ERC Advanced Grant PBDR-339251

if __name__ == "__main__":

  ## Overview:
  ExplanatoryMessage = '''This script removes the gap character "-" from sequences
  in a fasta file. Output is printed to stdout.'''
  
  
  
  # Define a function to check files exist, as a type for the argparse.
  def File(MyFile):
    if not os.path.isfile(MyFile):
      raise argparse.ArgumentTypeError(MyFile+' does not exist or is not a file.')
    return MyFile
  
  # Set up the arguments for this script
  ExplanatoryMessage = ExplanatoryMessage.replace('\n', ' ').replace('  ', ' ')
  parser = argparse.ArgumentParser(description=ExplanatoryMessage)
  parser.add_argument('FastaFile', type=File)
  parser.add_argument('-?', '--q-mark', action='store_true', help="Remove the "\
  " '?' character (used in shiver consensuses to mean missing coverage) too.")
  parser.add_argument('-TE', '--trim-missing-ends', action='store_true', help='''
  Trims any "?", "N" or "n" characters from the ends of each sequence.''')
  parser.add_argument('-1', '--first-seq-only', action='store_true', \
  help='''Return only the first sequence from the fasta file.''')
  args = parser.parse_args()
  
  UngappedSeqs = []
  for seq in SeqIO.parse(open(args.FastaFile),'fasta'):
    seq.seq = ungap(seq.seq)
    if args.q_mark:
      seq.seq = ungap(seq.seq, "?")
    if args.trim_missing_ends:
      SeqAsStr = str(seq.seq)
      SeqAsStr = sub("^[Nn?]+", "", SeqAsStr)
      SeqAsStr = sub("[Nn?]+$", "", SeqAsStr)
      seq.seq = Seq.Seq(SeqAsStr)
    UngappedSeqs.append(seq)
    if args.first_seq_only:
      break
  
  if UngappedSeqs == []:
    print('No sequences found in', args.FastaFile+'. Quitting.', file=sys.stderr)
    exit(1)
  
  SeqIO.write(UngappedSeqs, sys.stdout, "fasta")
  
