#!/usr/bin/env python

"""
motifMap.py => Look at the localization of a motif

To-do:
	* Motif mapping to genome
	* GTF feature selection
	* bedtools intersection
	* Genic vs intergenic
	* Coding vs non-coding (intron, exon, UTRs)
	* Proximity of motif to 5' and 3' splice sites
	* GO enrichment
"""

from __future__ import division
from natsort import natsorted
from fileParser import yield_fasta, parse_GTF
from sequence import Sequence

import os
import re
import argparse
import subprocess

# bedtools intersectBed -a <GTF> -b <MOTIFS.BED> -wo

def map2genome(genome_fasta, motif, output_dir):

	motif = motif[0]
	
	# Regex for overlapping motif occurences
	p = re.compile('(?=({}))'.format(motif))

	hits = set()
	for record in yield_fasta(genome_fasta):

		# Forward search
		for m in p.finditer(record.seq):
			hits.add((record.id, m.start(), m.start()+len(motif), '+'))

		# Reverse search
		seq = Sequence(record.seq).get_reverse_complement()
		seq_len = len(record.seq)
		for m in p.finditer(seq):
			hits.add((record.id, seq_len-(m.start()+len(motif)), seq_len-m.start(), '-'))

	with open('{}/{}.bed'.format(output_dir, motif), 'w') as fout: 
		for hit in natsorted(hits):
			c, s, e, strand = hit
			fout.write( '{}\t{}\t{}\t{}\t1000\t{}\n'.format(c, s, e, motif, strand) )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-d', '--dir', required=True, help="Output directory.")
	subparsers = parser.add_subparsers(dest='command', help="Sub-command help.")

	parser_a = subparsers.add_parser('map2genome', help="Map motif(s) to the genomic fasta sequence.")
	parser_a.add_argument('--genome', required=True, help="Genome fasta sequence.")
	parser_a.add_argument('--motif', nargs='+', required=True, help="List of motifs.")

	#parser_b.add_argument('-g', '--gtf', required=True, help="Transcript GTF annotation.")
	args = parser.parse_args()

	if args.command == 'map2genome':
		map2genome(args.genome, args.motif, args.dir)