#!/usr/bin/env python

"""
Get all alternative splicing events from a GTF file.
Replaces the get_AS_landscape.py script.
"""

from __future__ import print_function
import os
import sys
import argparse
from natsort import natsorted
from fileParser import parse_GTF

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

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

def get_SS(transcripts):
	
	import itertools

	# Throw all exons on one big heap and
	# ignore first and last exons (FOR NOW...)
	alt = { '{}-{}'.format(s, e): [[ s, e ]] for s, e in list(itertools.chain(*[ natsorted(transcripts[t_id]['exons'])[1:-1] for t_id in transcripts ])) }

	# Save exon and transcript relation
	c, strand = [ transcripts[transcripts.keys()[0]][x] for x in ['chr', 'strand'] ]
	exonmap = {}
	for t_id in transcripts:
		for s, e in transcripts[t_id]['exons']:
			try: exonmap['{}-{}'.format(s, e)].add(t_id)
			except KeyError: exonmap['{}-{}'.format(s, e)] = set([ t_id ])

	# Look for overlaps in the remaining exons
	for exon in natsorted(alt):
		s, e = map(int, exon.split('-'))

		for other_exon in alt:
			if exon != other_exon:
				x, y = map(int, other_exon.split('-'))
				if any([ s <= x <= e and e <= y, x <= s and s <= y <= e ]): alt[exon].append([ x, y ])

	# If there are any overlaps, look at those
	done = set()
	ir_exons = []
	for iso in natsorted(alt):

		# Remove intron retention exons
		if len(alt[iso]) > 2:

			# Get all the overlapping exons
			ir_check = {}
			for s, e in alt[iso]:
				for x, y in alt[iso]:
					if (s, e) != (x, y) and is_overlapping(s, e, x, y):
						try: ir_check['{}-{}'.format(s, e)].append([x, y])
						except KeyError: ir_check['{}-{}'.format(s, e)] = [[x, y]]

			''' 
			    ---XXXXXXXXXXXXXXXXXXXXX---

				---YYYY-------------ZZZZ---

				Check for overlapping exons. For the example above,
				remove exon 'X' if the exons it overlaps with, 'Y' 
				and 'Z', do not overlap with each other. 
			'''

			for exon in ir_check:
				if len(ir_check[exon]) > 1:
					s, e = ir_check[exon][0]
					for x, y in ir_check[exon][1:]:
						if not is_overlapping(s, e, x, y):
							ir_exons.append(map(int, exon.split('-')))

		if len(alt[iso]) > 1:
			for i in xrange(1, len(alt[iso])):

				try: 

					# Only analyse "simple" events where one of the splice sites is shared
					shrd = set( alt[iso][0] ).intersection( set( alt[iso][i] ) ).pop() # Shared coordinate
					iso1, iso2 = set( alt[iso][0] ).symmetric_difference( set( alt[iso][i] ) ) # Get the alternative coordinates
					index = alt[iso][0].index(shrd)

					if all([ alt[iso][0] not in ir_exons, alt[iso][i] not in ir_exons ]):

						# Check for duplicates
						if frozenset(sorted([shrd, iso1, iso2])) not in done:
							done.add(frozenset(sorted([shrd, iso1, iso2])))

							if index:
								ss, one, two = 'A3', '{}-{}'.format(iso1, shrd), '{}-{}'.format(iso2, shrd)
							else:
								ss, one, two = 'A5', '{}-{}'.format(shrd, iso1), '{}-{}'.format(shrd, iso2)

							event_id = '{}:{}:{}:{}:{}'.format(ss, c, one, two, strand)
							print( '{0}\tgtf2AS\texon\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}:alternative1"; transcript_id "{4}:alternative1"; transcripts "{5}";'.format(c, one.split('-')[0], one.split('-')[1], strand, event_id, ','.join(exonmap[one])) )
							print( '{0}\tgtf2AS\texon\t{1}\t{2}\t.\t{3}\t.\tgene_id "{4}:alternative2"; transcript_id "{4}:alternative2"; transcripts "{5}";'.format(c, two.split('-')[0], two.split('-')[1], strand, event_id, ','.join(exonmap[two])) )

				except KeyError: pass # No overlaps

def in_bounds(s, e, exons):

	left = min([ min(x, y) for x, y in exons ])
	right = max([ max(x, y) for x, y in exons ])
	return is_overlapping(s, e, left, right)

def get_ES(transcripts):
	
	first_id = transcripts.keys()[0]
	c, strand = [ transcripts[first_id][x] for x in ['chr', 'strand'] ]
	es_events = {}

	for t_id in natsorted(transcripts):
		for s, e in transcripts[t_id]['exons'][1:-1]:
			es = '{}-{}'.format(s, e)
			for other_t in transcripts:
				if t_id != other_t:

					# If an exon is not overlapping with any exon from another isoform and
					# this particular exon is within bounds of that isoform, it must
					# have been skipped!
					if not any([ is_overlapping(s, e, x, y) for x, y in transcripts[other_t]['exons'] ]) and in_bounds(s, e, transcripts[other_t]['exons']):
						try: es_events[es].add(other_t)
						except KeyError: es_events[es] = set([other_t])

	for es in natsorted(es_events):
		print('{0}\tgtf2AS\texon\t{1}\t{2}\t.\t{3}\t.\tgene_id "ES:{0}:{1}-{2}:{3}"; transcript_id "ES:{0}:{1}-{2}:{3}"; transcripts "{4}";'.format(c, es.split('-')[0], es.split('-')[1], strand, ','.join(natsorted(es_events[es]))))

def get_IR(transcripts):
	
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
					matches = [ [x, y] for x, y in transcripts[other_t]['exons']  if is_overlapping(s, e, x, y) ]
					if len(matches) > 1:

						# Check whether the overlapping exons are actually not 
						# in fact alternative 3' or 5' spliced exons
						a, b = matches[0]
						if any([ not is_overlapping(a, b, c, d) for c, d in matches[1:] ]):

							# Only select intron retention events with the same boundaries
							left = min([ min(x, y) for x, y in matches ])
							right = max([ max(x, y) for x, y in matches ])
							if natsorted([left, right]) == natsorted([s, e]):

								try: ir_events[event_id].add(t_id)
								except KeyError: ir_events[event_id] = set([ t_id ])

	for ir in natsorted(ir_events):
		print( '{0}\tgtf2AS\texon\t{1}\t{2}\t.\t{3}\t.\tgene_id "IR:{0}:{1}-{2}:{3}"; transcript_id "IR:{0}:{1}-{2}:{3}"; transcripts "{4}";'.format(chromosome, ir.split('-')[0], ir.split('-')[1], strand, ','.join(natsorted(ir_events[ir]))) )

def run(gtf_file, event_type):

	gtf = parse_GTF(gtf_file, select_feature="exon", get_introns=False)

	# Group per gene
	groups = {}
	for t in gtf:
		g_id = t.split('.')[0]
		try: groups[g_id].append(t)
		except KeyError: groups[g_id] = [t]

	for g_id in natsorted(groups):

		# Check if the gene has multiple transcripts
		if len(groups[g_id]) > 1:

			transcripts = { t_id: gtf[t_id] for t_id in groups[g_id] }
			if 'SS' in event_type: get_SS(transcripts)
			if 'ES' in event_type: get_ES(transcripts)
			if 'IR' in event_type: get_IR(transcripts)

if __name__ == '__main__':

	AS_choices = ['ES', 'SS', 'IR']

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Transcripts in GTF format.")
	parser.add_argument('-e', '--event-type', nargs='+', choices=AS_choices, default=AS_choices, help="AS events to look for.")
	parser.add_argument('--version', action='version', version='v0.3.0')
	args = parser.parse_args()

	run(args.gtf, args.event_type)