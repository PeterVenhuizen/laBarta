#!/usr/bin/env python

"""
Split fastQ reads in k-sized k-mers.
"""

from multiprocessing import Process, Queue
import gzip
import argparse
from itertools import izip_longest

import sys
sys.path.append('/home/venhuip8/scripts')
from fileParser import parseFastQ

def split_reads(queue, records, k=50):
	
	sub_reads = []
	for r in records:
		if r == None: break
		else:

			r_len = len(r[1])

			# If read can be split up exactly in n*k parts
			if (r_len % k) == 0:
				n_parts = r_len / k
				for i in xrange(n_parts): sub_reads.append( '>{}|{}\n{}\n'.format( r[0], i+1, r[1][(i*k):(i+1)*k] ) )

			# If read is too small to split in two non-overlapping parts
			elif ( r_len < (k*2) ) and r_len > k:
				sub_reads.append( '>{}|a\n{}\n'.format( r[0], r[1][:k] ) )
				sub_reads.append( '>{}|b\n{}\n'.format( r[0], r[1][(r_len-k):] ) )

			# If read is bigger than 2*k, but cannot be split-up in k-sized non-overlapping parts
			else:

				# Split into as many as possible k-size, non-overlapping parts
				n = (r_len // k) - 1
				for i in xrange(n):
					sub_reads.append( '>{}|{}\n{}\n'.format( r[0], i+1, r[1][(i*k):(i+1)*k] ) )

				# Get the last two overlapping parts
				part = r[1][(i+1)*k:]
				sub_reads.append( '>{}|a\n{}\n'.format( r[0], part[:k] ) )
				sub_reads.append( '>{}|b\n{}\n'.format( r[0], part[(len(part)-k):] ) )

	queue.put(''.join(sub_reads))

def grouper(fastQIter, n, fillvalue=None):
	""" Collect data into fixed-length chunks of blocks.
		grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx """

	args = [fastQIter] * n
	return izip_longest(fillvalue=fillvalue, *args)

def get_chunks(lst, n):
	return [ lst[i::n] for i in xrange(n) ]

if __name__ == '__main__':

	# Parse arguments
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-i', '--input', help="input fastQ file", required=True)
	parser.add_argument('-o', '--output', help="Output fasta file", required=True)
	parser.add_argument('-k', '--kmer-size', help="K-mer size (default = 50)", default=50, nargs='?', type=int)
	parser.add_argument('-p', '--processes', help="Number of parallel processes (default = 4)", default=4, nargs='?', type=int)
	args = parser.parse_args()
	
	# Get fastQ iterator
	fastQIter = parseFastQ(args.input)

	# Open output stream
	fout = gzip.open(args.output, 'wt') if args.output.endswith('gz') else open(args.output, 'w')

	queue = Queue()
	for records in grouper(fastQIter, args.processes * 16384):

		chunks = get_chunks(records, args.processes)

		processes = [ Process(target=split_reads, args=(queue, x)) for x in chunks ]

		for p in processes:
			p.start()

		#for p in processes:
		#	p.join()

		results = [ queue.get() for p in processes ]

		fout.write(''.join(results))

	fout.close()