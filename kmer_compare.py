#!/usr/bin/env python

"""
kmer_compare.py -> Compare the kmer frequency/enrichment
between two groups of sequences. Compare to frequency to
expected frequency based on the base composition and calculate
enrichment using the Fisher's Exact test.
"""

from __future__ import division
from fileParser import yield_fasta
from sequence import Sequence

from collections import Counter
from natsort import natsorted
from scipy.stats import fisher_exact
import argparse

def prep_analysis(fasta_file, kmers):

	ACGT, kmers_freq, sizes = Counter(), {}, []
	for record in yield_fasta(fasta_file):
		for b in record.seq:
			ACGT[b] += 1

		sizes.append(len(record.seq))

	# Count number of hits/misses per kmer
	for record in yield_fasta(fasta_file):
		for k in kmers:
			i = str(len(k))
			
			if i not in kmers_freq: kmers_freq[i] = {}
			if k not in kmers_freq[i]: kmers_freq[i][k] = { 'hit': 0, 'miss': 0 }

			if k in record.seq: kmers_freq[i][k]['hit'] += 1
			else: kmers_freq[i][k]['miss'] += 1

	return { b: ACGT[b]/sum(ACGT.values()) for b in ACGT }, kmers_freq, sizes

def kmer_compare(sample_fasta, background_fasta):

	koi = set() # kmers of interest
	for record in yield_fasta(sample_fasta):
		#for k in xrange(6, 21):
		for k in [5]:
			N = (len(record.seq)-k)+1
			for i in xrange(0, N, 1):
				koi.add(record.seq[i:(i+k)])

	#koi = set(['ACAG', 'AGAG', 'CAGA', 'CAG', 'GACAGA', 'ACAGAGA', 'ATTCAGA', 'ACAGAAG', 'AGACAGA', 'GACAGAA', 'CAGACAA', 'ATCAGAC', 'CAGCAGA', 'ACAGACA', 'CAACCAG', 'CAGACGA', 'AGCAGAC', 'CAGACAG', 'ACAGACG', 'GACAGAC'])

	print "I've got my kmers :)"

	smpl_freq, smpl_kmers, smpl_sizes = prep_analysis(sample_fasta, koi)
	print 'looked at my sample'
	bkgd_freq, bkgd_kmers, bkgd_sizes = prep_analysis(background_fasta, koi)
	print 'and checked the background'

	for i in natsorted(smpl_kmers):
		for k in natsorted(smpl_kmers[i]):

			smpl_hit, smpl_miss = [ smpl_kmers[i][k][x] for x in ['hit', 'miss'] ]
			bkgd_hit, bkgd_miss = [ bkgd_kmers[i][k][x] for x in ['hit', 'miss'] ]

			oddsratio, pvalue = fisher_exact([[smpl_hit, bkgd_hit], [smpl_miss, bkgd_miss]], alternative='greater')
			#if pvalue < 0.05 and smpl_hit > (len(smpl_sizes)/10):
			if pvalue < 0.05:
				print '{}\t{}\t{}\t{}\t{}\t{:.7f}'.format(k, smpl_hit, smpl_miss, bkgd_hit, bkgd_miss, pvalue)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-s', '--sample', required=True, help="Sample fasta file.")
	parser.add_argument('-b', '--background', required=True, help="Background fasta file.")
	args = parser.parse_args()

	kmer_compare(args.sample, args.background)