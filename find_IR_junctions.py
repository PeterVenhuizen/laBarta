#!/usr/bin/env python

"""
Look for potential intron retention events. An IR events is considered if:
	1) 	All bases of the junction are covered by reads
	2) 	The median coverage of the junction must be above min_cov (default = 5) and
		be at least 5% of the junction supporting reads

	Created on: 17-09-2015 (Peter Venhuizen)
	Last edited: 22-02-2017 (Peter Venhuizen)
"""

from __future__ import division
from fileParser import yield_junctions
from natsort import natsorted
import numpy as np
import argparse
import subprocess

def find_IR_junctions(junctions_bed, aligned_reads_bam, genome_fa_index, min_cov=5):
	
	# Parse junctions and create a junction index
	junctions, junc_indx = {}, {}
	for j in yield_junctions(junctions_bed):

		c, s, e, depth, strand = [ j[x] for x in ['chr', 'junc_start', 'junc_end', 'depth', 'strand'] ]
		locus = '{}:{}-{}'.format(c, s, e)
		junctions[locus] = { 'chr': c, 'start': s, 'end': e, 'strand': strand, 'depth': depth, 'length': (e-s)+1 }

		if c not in junc_indx: junc_indx[c] = {}
		if strand not in junc_indx[c]: junc_indx[c][strand] = {}
		for i in xrange(s, e+1):
			try: junc_indx[c][strand][str(i)].append(locus)
			except KeyError: junc_indx[c][strand][str(i)] = [locus]

	# Get coverage
	IR_junc = {}
	for strand in ['-', '+']:

		# Generate genome coverage bed
		cmd = "genomeCoverageBed -ibam {} -g {} -d -split -strand {} > tmp.genomeCoverageBed".format(aligned_reads_bam, genome_fa_index, strand)
		subprocess.call(cmd, shell=True)

		# Allocate coverage
		junc_cov = {}
		for line in open("tmp.genomeCoverage.bed"):
			c, pos, cov = line.rstrip().split('\t')
			if c in junc_indx:
				if pos in junc_indx[c][strand]:
					for locus in junc_indx[c][strand][pos]:
						try: junc_cov[locus].append( float(cov) )
						except KeyError: junc_cov[locus] = [ float(cov) ]
		subprocess.call( "rm tmp.genomeCoverage.bed", shell=True )

		# Filter junctions
		for locus in junc_cov:

			min_req_cov = junctions[locus]['depth'] * 0.05 # 5% of junction supporting reads
			median_cov = np.median(junc_cov[locus])

			if min(junc_cov[locus]) > 0 and median_cov > min_cov and median_cov > min_req_cov: IR_junc[locus] = median_cov

	# Output
	print '#chr\tstart\tend\tstrand\tsize\tmod3\tmedian_cov'
	for l in natsorted(IR_junc):
		print '{}\t{}\t{}\t{}\t{}\t{}'.format(junctions[l]['chr'], junctions[l]['start'], junctions[l]['end'], junctions[l]['length'], junctions[l]['length'] % 3, IR_junc[l])

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-j', '--junctions', required=True, help="TopHat2 junctions.bed file.")
	parser.add_argument('-b', '--bam', required=True, help="TopHat2 aligned reads bam file.")
	parser.add_argument('--genome-index', default="/home/venhuip8/data/master_transcriptome/a_thaliana/a_thaliana_TAIR10.fa.fai", help="Genome fasta index for genomeCoverageBed")
	parser.add_argument('--min-coverage', default=5, type=int, help="Minimum intron coverage")
	args = parser.parse_args()

	find_IR_junctions(args.junctions, args.bam, args.genome_index, args.min_coverage)