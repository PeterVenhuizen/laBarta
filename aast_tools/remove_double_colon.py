#!/usr/bin/env python

import argparse

def yield_fasta(f):
	''' Simple fasta parser that yield's the identifier and sequence of each record '''
	
	class SeqRecord:
		def __init__(self, seq_id, seq):
			self.id = seq_id
			self.seq = seq

	seq = ''
	for line in f:
		if line.startswith('>'):
			if len(seq):
				yield(SeqRecord(identifier, seq))
				seq = ''
			identifier = line[1:].rstrip()
		else: seq += line.rstrip()
	yield(SeqRecord(identifier, seq))

def remove_double_colon(fasta_file):
	
	w = 80
	for record in yield_fasta(fasta_file):
		print('>{}'.format(record.id.split('::')[0]))
		for i in xrange(0, len(record.seq), w):
			print(record.seq[i:i+w])

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--fasta', type=argparse.FileType('r'), default='-', help="Exon skipping events.")
	args = parser.parse_args()

	remove_double_colon(args.fasta)