#!/usr/bin/env python

'''
	Name: get_AS_landscape.py
	Author: Petrus Andreas Venhuizen
	E-mail: peter.venhuizen@univie.ac.at
	Institute: Max F. Perutz Laboratories / Medical University of Vienna
	Description: Report the alternative splicing landscape of a set of transcripts.
	Required packages: natsort (https://pypi.python.org/pypi/natsort)
	
	Created on: 01-10-2015
	Last changed: 22-02-2017
	
	To-do:
	- 	Change parse_gtf(). Allow for missing exon_number attribute.
	-	Add gff file support
'''

from __future__ import print_function
import os
import argparse
from natsort import natsorted
from fileParser import parse_GTF

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def get_longest_ref(trs):
	""" Select the longest transcript """

	import operator
	l = { t_id: sum([ abs(e-s) for s, e in trs[t_id]['exons'] ]) for t_id in trs }
	return max(l.iteritems(), key=operator.itemgetter(1))[0]

def find_unknown_start(start, end, trs):
	if start > end: return sum([ 1 for s, e in trs if start < e ])
	else: return sum([ 1 for s, e in trs if end > s ])

def match_isoform_exons(trs, ref, strand):
	
	matches, matched_ref, ref_c = [], set([]), {}
	i = 0
	for s, e in trs:
		
		if strand == '-': s, e = e, s
		exon, ref_exons = [s, e], []
		if s > e: s, e = e, s

		for S, E in ref:
			if strand == '-': S, E = E, S
			if S > E: x, y = E, S
			else: x, y = S, E

			add = False
			if s <= x <= e and s <= y <= e: add = True # Ref exon retained in trs exon
			elif x <= s <= y and x <= e <= y: add = True # Trs exon retained in ref exon
			elif x <= s <= y and x < e: add = True # Right overhang overlap with ref exon
			elif s < x and x <= e <= y: add = True # Left overhang overlap with ref exon

			if add:
				ref_exons.append([S, E])
				matched_ref.add((S, E))
				try: ref_c['{}-{}'.format(S, E)].append(i)
				except KeyError: ref_c['{}-{}'.format(S, E)] = [i]

		i += 1
		matches.append([exon, ref_exons])

	# Check if any of the reference exons are not assigned
	if len(ref) != len(matched_ref):
		leftovers = [ [s, e] for s, e in ref if (s, e) not in matched_ref ]
		for s, e in leftovers:
			i = find_unknown_start(s, e, trs)

			# In case of more skipped exons in the transcript
			i += sum([ 1 for j in xrange(len(matches)) if not sum(matches[j][0]) ])

			try: matches[i:i] = [[[0, 0], [[s, e]]]]
			except IndexError: matches.append([[[0, 0], [[s, e]]]])

	return matches, ref_c

def go_landscaping(transcripts, ref_id, event_type, detailed=False):
	
	# Get reference information
	ref_t = transcripts.pop(ref_id)
	n_ref_exons = len(ref_t['exons'])

	# Analyze each transcript
	t_list = transcripts.keys()
	for t_id in natsorted(t_list):

		# Required initialization variables
		i = 0
		n_trs_exons = len(transcripts[t_id]['exons'])
		alt_tss, alt_tts, alt_intron = False, False, False
		ref_count = 0
		landscape = {}
		ref_exon_skip = {}
		ref_IR = {}
		ref_SS = {}
		bed_track = []
		print_list = [[], [], []]

		# Try to match the isoform exons to the reference transcript
		matches, ref_c = match_isoform_exons(transcripts[t_id]['exons'], ref_t['exons'], transcripts[t_id]['strand'])
		for exon, ref_exons in matches:

			s, e = exon # Get exon start and end
			trs_has_exon = sum([s, e]) # Check if it is exon skipping
			n_matching_ref = len(ref_exons) # Number of matching ref exons for this trs exon
			bed_start, bed_end = s, e # For later bed-track output

			#print(i, ref_count, exon, ref_exons, s, e, bed_start, bed_end)

			if trs_has_exon or n_matching_ref and i < n_trs_exons:

				if n_matching_ref: ref_count += n_matching_ref # For exon skipping

				if 'ES' in event_type:
					# 1) Cassette exon / exon skipping
					# Skipping of reference exon
					if not trs_has_exon:

						if i and i < n_trs_exons and ref_count > 1:
							diff = abs(ref_exons[0][0]-ref_exons[0][1])+1
							bed_start, bed_end = ref_exons[0]
							ref_exon_skip[str(ref_count)] = 'ALT_EXON(-{})'.format(diff) if detailed else 'ALT_EXON'

					# Inclusion of non-reference exon
					elif not n_matching_ref and 1 <= i < n_ref_exons and ref_count >= 1:
						if i+1 < n_trs_exons: # If it is not the last exon in the transcript

							diff = abs(s-e)+1
							try: landscape[str(i+1)].append( 'ALT_EXON(+{})'.format(diff) if detailed else 'ALT_EXON' )
							except KeyError: landscape[str(i+1)] = [ 'ALT_EXON(+{})'.format(diff) if detailed else 'ALT_EXON' ]

				if 'FL' in event_type:
					# 2) Alternative transcription start site (ALT_TSS), first exon
					# starts before or after the first reference transcript exon
					if not i:

						# Check if an reference exon was found
						if not alt_tss:
							if (n_matching_ref and s != ref_exons[0][0]) or not n_matching_ref:
								alt_tss = True
								try: landscape['1'].append('ALT_TSS')
								except KeyError: landscape['1'] = ['ALT_TSS']

				if 'SS' in event_type:
					# 3) Alternative 3' splice site (ALT_3SS)
					if i and ref_count > 1: # Skip the first exon

						if n_matching_ref and trs_has_exon and not alt_intron and not alt_tss:
							if s != ref_exons[0][0]:

								# Calculate transcript size difference
								diff = ref_exons[0][0] -s if s < e else s - ref_exons[0][0]
								diff = '+{}'.format(diff) if diff > 0 else str(diff)

								try: landscape[str(i+1)].append( 'ALT_3SS({})'.format(diff) if detailed else 'ALT_3SS' )
								except KeyError: landscape[(str(i+1))] = [ 'ALT_3SS({})'.format(diff) if detailed else 'ALT_3SS' ]

								ref_SS[str(ref_count)] = 'ALT_3SS'

				if 'IR' in event_type:
					# 4a) Reference intron retention / Potential exitron
					if n_matching_ref and trs_has_exon:

						# Check if multiple trs exons map to the ref exon
						r = '{}-{}'.format(*ref_exons[0])
						if len(ref_c[r]) > 1 and i:
							bed_end = transcripts[t_id]['exons'][ref_c[r][-1]][1]
							if not ref_c[r].index(i):

								if transcripts[t_id]['strand'] == '+':
									bed_start = transcripts[t_id]['exons'][ref_c[r][0]][1]
									bed_end = transcripts[t_id]['exons'][ref_c[r][-1]][0]
								else:
									bed_start = transcripts[t_id]['exons'][ref_c[r][0]][0]
									bed_end = transcripts[t_id]['exons'][ref_c[r][-1]][1]

								# Calculate nt lost compared to the ref
								diff = 0
								for j in xrange(1, len(ref_c[r])):
									k = ref_c[r][j]
									s1, e1 = transcripts[t_id]['exons'][k-1]
									s2, e2 = transcripts[t_id]['exons'][k]
									diff += abs(e1-s2)

								ref_IR[str(i+1)] = 'POT_EI_{}-{}({})'.format(ref_c[r][0]+1, ref_c[r][-1]+1, diff) if detailed else 'POT_EI_{}-{}'.format(ref_c[r][0]+1, ref_c[r][-1]+1)
								alt_intron = True

					# 4b) Intron retention
					if n_matching_ref > 1:

						# Calculate transcript length difference
						diff = 0
						for j in xrange(1, len(ref_exons)):
							s1, e1 = ref_exons[j-1]
							s2, e2 = ref_exons[j]
							diff += abs(e1-s2)

						try: landscape[str(i+1)].append('IR_{}-{}(+{})'.format(i+1, i+n_matching_ref, diff) if detailed else 'IR_{}-{}'.format(i+1, i+n_matching_ref))
						except KeyError: landscape[str(i+1)] = [ 'IR_{}-{}(+{})'.format(i+1, i+n_matching_ref, diff) if detailed else 'IR_{}-{}'.format(i+1, i+n_matching_ref) ]

				if 'SS' in event_type:
					# 5) Alternative 5' splice site (ALT_5SS)
					if i+1 != n_trs_exons and not alt_intron and not alt_tss: # Skip the last exons and alternative introns
						
						if n_matching_ref and trs_has_exon:
							if e != ref_exons[-1][1]:

								# Calculate transcript length difference
								diff = e - ref_exons[-1][1] if s < e else ref_exons[-1][1] - e
								diff = '+{}'.format(diff) if diff > 0 else str(diff)

								try: landscape[str(i+1)].append('ALT_5SS({})'.format(diff) if detailed else 'ALT_5SS')
								except KeyError: landscape[str(i+1)] = [ 'ALT_5SS({})'.format(diff) if detailed else 'ALT_5SS' ]

								#ref_SS[str(ref_count)] = 'ALT_5SS'

				if 'FL' in event_type:
					# 6) Alternative transcription termination site (ALT_TTS), last
					# exon ends before or after the last reference transcript exon
					if i+1 == n_trs_exons and trs_has_exon:

						# Check if an reference exon was found
						if not alt_tts:

							# If it overlaps with an reference exon or not
							if (n_matching_ref and e != ref_exons[-1][1]) or not n_matching_ref: 
								alt_tss = True
								try: landscape[str(i+1)].append('ALT_TTS')
								except KeyError: landscape[str(i+1)] = ['ALT_TTS']

				if trs_has_exon: i += 1

				# Build bed-track output lines
				c, strand = [ transcripts[t_id][x] for x in ['chr', 'strand'] ]
				if bed_start > bed_end: bed_start, bed_end = bed_end, bed_start

				if str(i) in landscape and trs_has_exon and i not in print_list[0]:
					print_list[0].append(i)
					bed_track.append( '{}\t{}\t{}\t{}|{}:{}\t1000\t{}'.format(c, bed_start, bed_end, t_id, i, ','.join(landscape[str(i)]), strand) )

				if str(i) in ref_IR and i not in print_list[0]:
					print_list[0].append(i)
					bed_track.append( '{}\t{}\t{}\t{}|{}:{}\t1000\t{}'.format(c, bed_start, bed_end, ref_id, i, ref_IR[str(i)], strand) )

		yield('\n'.join(bed_track))

def run(gtf_file, ref, event_type, output_dir):
	
	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Parse transcripts in the gtf_file
	gtf = parse_GTF(gtf_file, select_feature="exon", get_introns=False)

	# Group per gene
	groups = {}
	for t in gtf:
		g_id = t.split('.')[0]
		try: groups[g_id].append(t)
		except KeyError: groups[g_id] = [t]

	# Parse reference transcripts
	reference = { line.rstrip().split('.')[0]: line.rstrip() for line in open(ref) }

	for g_id in natsorted(groups):

		# Check if the gene has multiple transcripts
		if len(groups[g_id]) > 1:

			transcripts = { t_id: gtf[t_id] for t_id in groups[g_id] }

			# Select the longest transcript as the reference if there
			# is not reference in the file
			if g_id not in reference:
				new_ref = get_longest_trs(transcripts)
				reference[g_id] = new_ref
				eprint("WARNING: Missing reference transcript for gene {}. Selected longest transcript {} as reference.".format(g_id, new_ref))

			# Check if the reference transcript is in the GTF file
			# and run if everything is fine
			try:
				for bed_output in go_landscaping(transcripts, reference[g_id], event_type, True):
					print(bed_output)

			except KeyError:
				eprint("Warning: Reference transcript {} not found in the GTF file. Skipping {} transcripts.".format(reference[g_id], len(transcripts)))

if __name__ == '__main__':

	AS_choices = ['ES', 'SS', 'IR', 'FL']

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Transcripts in GTF format.")
	parser.add_argument('-r', '--reference-trs', required=True, help="List of reference transcripts for the input GTF genes.")
	parser.add_argument('-e', '--event-type', nargs='+', choices=AS_choices, default=AS_choices, help="AS events to look for.")
	parser.add_argument('-o', '--output', default="./", help="Output path/file")
	parser.add_argument('--version', action='version', version='v0.3.0')
	args = parser.parse_args()

	run(args.gtf, args.reference_trs, args.event_type, args.output)