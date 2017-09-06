#!/usr/bin/env python

"""
mxi.py => Detect Mutually Exclusive Introns in a given GTF file. 

	---|||||-----|||||||||||||||---
	---|||||||||||||||-----|||||---
	---|||||-----|||||-----|||||---

	AStalavista MXI code => 1^2-,3^4-

"""

import argparse
import itertools
from fileParser import parse_GTF
from natsort import natsorted

def is_overlapping(s, e, x, y):
	''' Check if s-e overlaps with x-y '''

	s, e = natsorted((s, e))
	x, y = natsorted((x, y))

	return any([
		s <= x <= e and s <= y <= e, # Ref exon retained in trs exon
		x <= s <= y and x <= e <= y, # Trs exon retain in ref exon
		x <= s <= y and x < e, # Right overlap with ref exon
		s < x and x <= e <= y # Left overlap with ref exon
	])

def is_contained_in(s, e, x, y):
	''' Check if s-e is contained in x-y '''

	s, e = natsorted((s, e))
	x, y = natsorted((x, y))
	return x < s < y and x < e < y

def mxi(gtf_file):

	gtf = parse_GTF(gtf_file, select_feature="exon")

	# Group per gene
	groups = {}
	for t in gtf:
		g_id = t.split('.')[0]
		try: groups[g_id].append(t)
		except KeyError: groups[g_id] = [t]

	for g_id in natsorted(groups):

		if len(groups[g_id]) > 1:

			transcripts = { t_id: gtf[t_id] for t_id in groups[g_id] }
			first_id = transcripts.keys()[0]
			chromosome, strand = [ transcripts[first_id][x] for x in ['chr', 'strand'] ]
			ir_exons = {}

			# Get IR exons
			for t_id in natsorted(transcripts):
				for s, e in transcripts[t_id]['exons']:

					# Look for overlapping exons in other isoforms
					exon_id = '{}-{}'.format(s, e)
					for other_t in transcripts:
						if t_id != other_t:

							# Get all the overlapping exons
							matches = [ [x, y] for x, y in transcripts[other_t]['exons'] if is_overlapping(s, e, x, y) ]
							if len(matches) > 1:

								# Check whether the overlapping exons are actually not
								# in fact alternative 5' and 3' spliced exons
								a, b = matches[0]
								if any([ not is_overlapping(a, b, c, d) for c, d in matches[1:] ]):
									try: ir_exons[exon_id]['trs'].add(t_id)
									except KeyError: ir_exons[exon_id] = { 'trs': set([ t_id ]), 'introns': set() }
			ir_list = [ map(int, ir.split('-')) for ir in ir_exons ]

			# Assign the retained introns
			introns = {}
			for t_id in natsorted(transcripts):
				for s, e in transcripts[t_id]['introns']:

					# Save intron/transcript info
					intron_id = '{}-{}'.format(s, e)
					try: introns[intron_id].add(t_id)
					except KeyError: introns[intron_id] = set([ t_id ])

					for x, y in ir_list:
						if is_contained_in(s, e, x, y): ir_exons['{}-{}'.format(x, y)]['introns'].add((s, e))

			# Look for overlapping IR events
			while len(ir_list) >= 1:
				s, e = ir_list.pop()
				matches = [ [x, y] for x, y in ir_list if is_overlapping(s, e, x, y) ]

				if len(matches):

					ir1, ir2 = [ '{}-{}'.format(x, y) for x, y in natsorted([ [s, e], [matches[0][0], matches[0][1]] ]) ]

					# Transform sets to lists
					ir_exons[ir1]['introns'] = list(ir_exons[ir1]['introns'])
					ir_exons[ir2]['introns'] = list(ir_exons[ir2]['introns'])

					# Iterate over the potentially multiple retained introns
					for i in ir_exons[ir1]['introns']:

						for j in ir_exons[ir2]['introns']:

							# Check that the retained introns are not the same and
							# that they share a common "middle" exon
							if ir_exons[ir1]['introns'] != ir_exons[ir2]['introns'] and i[1]+1 == int(ir2.split('-')[0]):

								print '{0}\t{1}\t{1};MXI:{0}:{2}:{3}-{4}:{5}-{6}:{7}:{8}\t{9}\t{10}'.format(chromosome, g_id, ir1.split('-')[0], i[0], i[1], j[0], j[1], ir2.split('-')[1], strand, ','.join(ir_exons[ir1]['trs']), ','.join(ir_exons[ir1]['trs']|ir_exons[ir2]['trs']))

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Transcripts in GTF format.")
	args = parser.parse_args()

	mxi(args.gtf)