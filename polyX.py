#!/usr/bin/env python

"""
Find polyX stretches in a fasta file.
"""

from fileParser import yield_fasta
from natsort import natsorted
import argparse

def polyX(fasta_file, X='A', nX=2, max_gap=0, min_len=2):
	
	# Check that min_len is not smaller than nX
	min_len = min_len if min_len > nX else nX

	for record in yield_fasta(fasta_file):

		# Initialize variables
		seq = record.seq.upper()
		seq_len = len(seq)
		start, end, gap = None, None, 0
		result = []

		for i in xrange(seq_len):

			# Equals X, so continue
			if seq[i] == X:
				if type(start) is int: end = i
				else: start = i
				gap = 0

			# No consecutive X's found
			else:

				gap += 1
				if gap > max_gap:

					# If start and end are defined and stretch contains nX X's
					if type(start) is int and type(end) is int \
					and X * nX in seq[start:end+1] and len(seq[start:end+1]) >= min_len:

						result.append((start, seq[start:end+1]))

					# Reset variables
					start, end, gap = None, None, 0

		# Check at the end of the sequence
		if type(start) is int and type(end) is int \
		and X * nX in seq[start:] and len(seq[start:]) >= min_len:
			result.append((start, seq[start:]))

		# Output
		for start, stretch in result: print '{}\t{}\t{}'.format(record.id, start, stretch)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--fasta', required=True, help="Fasta file to search in.")
	parser.add_argument('-b', '--base', choices=['A', 'C', 'G', 'T'], help="Poly [base] stretch to look for.")
	parser.add_argument('-m', '--min-Xs', default=2, nargs='?', type=int, help="Minimum required consecutive Xs.")
	parser.add_argument('-g', '--max-gap', default=0, nargs='?', type=int, help="Maximum allowed gap size.")
	parser.add_argument('-l', '--min-len', default=2, nargs='?', type=int, help="Minimum required poly stretch length.")
	args = parser.parse_args()

	polyX(args.fasta, args.base, args.min_Xs, args.max_gap, args.min_len)