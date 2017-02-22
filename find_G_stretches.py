#!/usr/bin/env python

"""
Find G-stretches in a fasta file.
"""

from fileParser import yield_fasta
from natsort import natsorted
import argparse

def find_G_stretches(fasta_file, nG=2, max_gap=0, min_len=2):

	# Check that min_len is not smaller than nG
	min_len = min_len if min_len > nG else nG

	for record in yield_fasta(fasta_file):

		# Initialize variables
		seq = record.seq.upper()
		seq_len len(seq)
		start, end, gap = None, None, 0
		result = []

		for i in xrange(seq_len):

			# Is G, so continue
			if seq[i] == 'G':
				if type(start) is int: end = i
				else: start = i
				gap = 0

			# No consecutive G's found
			else:
				
				gap += 1
				if gap > max_gap:
					
					# If start and end are defined and stretch contains nG G's
					if type(start) is int and type(end) is int \
					and 'G' * nG in seq[start:end+1] and len(seq[start:end+1]) >= min_len:
						
						result.append(( start, seq[start:end+1] ))
					
					# Reset variables
					start = None
					end = None
					gap = 0
		
		# Check at the end of the sequence
		if type(start) is int and type(end) is int \
		and 'G' * nG in seq[start:-1] and len(seq[start:-1]) >= min_len:
			result.append(( start, seq[start:-1] ))

		# Output
		for start, stretch in result: fout.write( '{}\t{}\t{}'.format(record.id, start, stretch) )

if __name__ == '__main__': 

	parser = argparse.ArgumentParser(description="Look for G-stretches in a fasta file.")
	parser.add_argument('-f', '--fasta-file', help="Fasta file to search in", required=True)
	parser.add_argument('-m', '--min-Gs', help="Minimum required consecutive Gs", default=2, nargs='?', type=int)
	parser.add_argument('-G', '--max-gap', help="Maximum gap between G-stretches", default=0, nargs='?', type=int)
	parser.add_argument('-l', '--min-len', help="Minimum required stretch length", default=2, nargs='?', type=int)
	args = parser.parse_args()
	
	run(args.fasta_file, args.output, args.min_Gs, args.max_gap, args.min_len)