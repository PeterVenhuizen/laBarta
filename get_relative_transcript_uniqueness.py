#!/usr/bin/env python

"""
Calculate the relative transcript sequence uniqueness. 
"""

from __future__ import division
import argparse
from fileParser import parse_GTF
from natsort import natsorted

def get_relative_transcript_uniqueness(gtf_file, select_feature="exon", t_id_attr="transcript_id", attr_sep=' "'):
	
	# Parse gtf
	gtf = parse_GTF(gtf_file, select_feature, t_id_attr, attr_sep, False)

	# Group per gene
	groups = {}
	for t in gtf:
		g_id = t.split('.')[0]
		try: groups[g_id].append(t)
		except KeyError: groups[g_id] = [t]

	for g_id in natsorted(groups):

		if len(groups[g_id]) > 1:
			trs = natsorted(groups[g_id])

			# Get "gene" start and end
			gene_start = min([ min([ s for s, e in gtf[t]['exons'] ]) for t in trs ])
			gene_end = max([ max([ e for s, e in gtf[t]['exons'] ]) for t in trs ])
			gene_len = abs(gene_end-gene_start)

			# Build exon/intron strings
			trs_strings = []
			for t in trs:
				intron_exon = list('I' * gene_len)
				for s, e in gtf[t]['exons']: intron_exon[abs(s-gene_start):abs(e-gene_start)] = list('E' * abs(e-s))
				trs_strings.append(intron_exon)

			#for i, t in enumerate(trs):
			#	print t, len(trs_strings[i]), ''.join(trs_strings[i])

			# Unique position count array
			counts = [0] * len(trs)
			for i in xrange(0, gene_len):
				pos_str = ''.join([ ei_string[i] for ei_string in trs_strings ])
				if pos_str.count('E') == 1: counts[pos_str.index('E')] += 1

			# Output unique percentage per transcript
			for i, t_id in enumerate(trs):
				trs_len = trs_strings[i].count('E')
				print '{}\t{}\t{:2.6f}'.format(t_id, trs_len, (counts[i] / trs_len))

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Input GTF file")
	parser.add_argument('--select-feature', default="exon", help="Feature to select in the GTF file")
	parser.add_argument('--transcript-attr', default='transcript_id', help="The transcript identifier attribute name")
	parser.add_argument('--attr-sep', default=' "', help="GFF attribute and attribute value separator")
	args = parser.parse_args()

	get_relative_transcript_uniqueness(args.gtf, args.select_feature, args.transcript_attr, args.attr_sep)