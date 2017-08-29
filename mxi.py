#!/usr/bin/env python

"""
mxi.py => Detect Mutually Exclusive Introns in a given GTF file. 

	---|||||-----|||||||||||||||---

	---|||||||||||||||-----|||||---

	---|||||-----|||||-----|||||---

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
		#if g_id == "AT1G02120":

			transcripts = { t_id: gtf[t_id] for t_id in groups[g_id] }
			first_id = transcripts.keys()[0]
			chromosome, strand = [ transcripts[first_id][x] for x in ['chr', 'strand'] ]
			ir_events = {}

			for t_id in natsorted(transcripts):
				for s, e in transcripts[t_id]['exons']:

					event_id = '{}-{}'.format(s, e)

					# Look for overlapping exons in other isoforms
					for other_t in transcripts:
						if t_id != other_t:

							# Get all the overlapping exons
							matches = [ [x, y] for x, y in transcripts[other_t]['exons'] if is_overlapping(s, e, x, y) ]
							if len(matches) > 1:

								# Check whether the overlapping exons are actually not
								# in fact alternative 3' and 5' spliced exons
								a, b = matches[0]
								if any([ not is_overlapping(a, b, c, d) for c, d in matches[1:] ]):

									try: ir_events[event_id].add(t_id)
									except KeyError: ir_events[event_id] = set([ t_id ])

			# Retrieve the retained introns
			all_introns = set([ (x, y) for x, y in list(itertools.chain(*[ transcripts[t_id]['introns'] for t_id in transcripts ])) ])
			all_IR = [ map(int, ir.split('-')) for ir in ir_events ]
			ret_introns = {}
			for s, e in all_IR:
				introns = [ [x, y] for x, y in all_introns if is_contained_in(x, y, s, e) ]
				if len(introns) == 1: ret_introns['{}-{}'.format(s, e)] = introns[0]

			# Look for overlapping IR events
			all_IR = [ map(int, ir.split('-')) for ir in ret_introns ]
			while len(all_IR) > 1:
				s, e = all_IR.pop()
				matches = [ [x, y] for x, y in all_IR if is_overlapping(s, e, x, y) ]
				if len(matches):

					ir1, ir2 = [ '{}-{}'.format(x, y) for x, y in natsorted([[s, e], [matches[0][0], matches[0][1]]]) ]

					if ret_introns[ir1] != ret_introns[ir2]:

						try: 
							print g_id
							print ir1, ret_introns[ir1]
							print ir2, ret_introns[ir2]
							print '{}\t{}:{}-{}:{}-{}:{}'.format(g_id, ir1.split('-')[0], ret_introns[ir1][0], ret_introns[ir1][1], ret_introns[ir2][0], ret_introns[ir2][1], ir2.split('-')[1])
							raw_input()
						except IndexError: pass

			#for ir in natsorted(ir_events):
			#	print ir, ir_events[ir]

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Transcripts in GTF format.")
	args = parser.parse_args()

	mxi(args.gtf)