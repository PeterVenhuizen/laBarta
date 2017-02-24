#!/usr/bin/env python

"""
Calculate the PSI/PIR/PSU values. See individual functions for more details.
"""

from __future__ import division
from __future__ import print_function
from natsort import natsorted
import argparse
import sys
from utils import parse_fasta_file, is_overlapping

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def parse_eff_count(eff_count):
	return { line.split('\t')[0]: int(line.rstrip().split('\t')[1]) for line in open(eff_count) }

def crc(EEJ_count, EEJ_mappability, EEJ_id, max_mappability=35):


	eff = EEJ_mappability[EEJ_id] if EEJ_id in EEJ_mappability else 1
	if EEJ_id in EEJ_count:
		return EEJ_count[EEJ_id] * (max_mappability/eff)
	return 0

def QC(inc, exc, AS_type):

	if AS_type in ["ES", "IR"]:
		if exc >= 100 or inc[0] >= 100 and inc[1] >= 100: QS = "SOK"
		elif 100 > exc >= 20 or any([ 100 > inc[0] >= 20 and 100 > inc[1] >= 15, 100 > inc[0] >= 15, 100 > inc[1] >= 20 ]): QS = "OK"
		elif 20 > exc >= 15 or any([ 20 > inc[0] >= 15 and 15 > inc[1] >= 10, 15 > inc[0] >= 10 and 20 > inc[1] >= 15 ]): QS = "LOW"
		elif 15 > exc >= 10 or any([ 15 > inc[0] >= 10 and 10 > inc[1] >= 5, 10 > inc[0] >= 5 and 15 > inc[1] >= 10 ]): QS = "VLOW"
		else: QS = "N"

	elif AS_type == "ME":
		if exc >= 100 or sum(inc) >= 100: QS = "SOK"
		elif 100 > exc >= 20 or 100 > sum(inc) >= 20: QS = "OK"
		elif 20 > exc >= 15 or 20 > sum(inc) >= 15: QS = "LOW"
		elif 15 > exc >= 10 or 15 > sum(inc) >= 10: QS = "VLOW"
		else: QS = "N"

	elif AS_type == "SS":
		if sum(inc) >= 100: QS = "SOK"
		elif 100 > sum(inc) >= 40: QS = "OK"
		elif 40 > sum(inc) >= 20: QS = "LOW"
		elif 20 > sum(inc) >= 10: QS = "VLOW"
		else: QS = "N"

	return QS

def RIS(inc):

	try: 
		r = max(inc) / min(inc)
		if r < 2: QS = "OK"
		elif 5 > r >= 2: QS = "B1"
		elif r >= 5: QS = "B2"
	except ZeroDivisionError: QS = "Bn"

	return QS

def ES(fasta_files, eff, count):
	''' Determine all (multi) exon skipping events from the gene/transcript/exon coordinates,
	identify all inclusion and exclusion EEJs and then calculate the PSI. '''

	# Parse transcripts
	seqs = {}
	for f in fasta_files:
		seqs = parse_fasta_file(seqs, f)

	print('#GENE\tES\tCOORD\tPSI\tQUALITY_SCORES\tEVENT_SIZE')
	for g_id in natsorted(seqs):

		exons = seqs[g_id]
		isFwd = exons.pop('isFwd')
		INC, ES = set(), {}

		# Collect the event information
		for t in natsorted(exons):

			ekeys = natsorted(exons[t].keys()) if isFwd else natsorted(exons[t].keys(), reverse=True)

			for x in ekeys[1:-1]:

				# Upstream donor/acceptor
				i = ekeys.index(x)
				es_start = x.split('-')[0] if isFwd else x.split('-')[1]

				sub = ekeys[i:-1]
				for y in sub:

					# Get exon skipping event
					j = ekeys.index(y)
					es_end = y.split('-')[1] if isFwd else y.split('-')[0]
					event = '{}-{}'.format(es_start, es_end)
					if event not in ES: ES[event] = { 's': int(es_start), 'e': int(es_end), 'CiA': set(), 'ACj': set(), 'INT': set(), 'EXC': set() }

					# Get exclusion EEJ
					donor = ekeys[i-1].split('-')[1] if isFwd else ekeys[i-1].split('-')[0]
					acceptor = ekeys[j+1].split('-')[0] if isFwd else ekeys[j+1].split('-')[1]
					ES[event]['EXC'].add('{}-{}'.format(donor, acceptor))

					# Get CiA and ACj
					CiA = '{}-{}'.format(donor, es_start)
					ACj = '{}-{}'.format(es_end, acceptor)
					ES[event]['CiA'].add(CiA)
					ES[event]['ACj'].add(ACj)

					# Save all inclusion EEJs. Needed to later add missing exclusion EEJs
					INC.update([CiA, ACj])

					# Get internal inclusion EEJ
					for k in xrange(i, j):
						e = ekeys[k].split('-')[1] if isFwd else ekeys[k].split('-')[0]
						s = ekeys[k+1].split('-')[0] if isFwd else ekeys[k+1].split('-')[1]
						ES[event]['INT'].add('{}-{}'.format(e, s))

		# Remove Alt3 and Alt5 events from ES
		# Remove any events which are overlapping
		#rm = set()
		#for event in ES:
		#	s, e = natsorted(map(int, event.split('-')))
		#	for other_event in ES:
		#		if other_event not in rm:
		#			x, y = natsorted(map(int, other_event.split('-')))
		#			if (s, e) != (x, y) and is_overlapping(s, e, x, y): rm.add(other_event)
		#for e in rm: del ES[e]

		# Add potentially missing exclusion EEJs and calculate PSI
		for event in natsorted(ES):

			es_start, es_end = ES[event]['s'], ES[event]['e']

			# Look for C1Cj
			for c in ES[event]['CiA']:
				c_start = c.split('-')[0]
				for d in INC:
					d_start, d_end = d.split('-')
					if c_start == d_start:
						if int(d_end) > es_end and es_start < es_end or int(d_end) < es_end and es_start > es_end: ES[event]['EXC'].add(d)

			# Look for CiC2
			for c in ES[event]['ACj']:
				c_end = c.split('-')[1]
				for d in INC:
					d_start, d_end = d.split('-')
					if c_end == d_end:
						if int(d_start) < es_start and es_start < es_end or int(d_start) > es_start and es_start > es_end: ES[event]['EXC'].add(d)

			# Get read counts
			CiA_raw = sum([ count[g_id+'_'+i] for i in ES[event]['CiA'] if g_id+'_'+i in count ])
			ACj_raw = sum([ count[g_id+'_'+i] for i in ES[event]['ACj'] if g_id+'_'+i in count ])
			EXC_raw = sum([ count[g_id+'_'+i] for i in ES[event]['EXC'] if g_id+'_'+i in count ])

			# Get corrected read counts
			CiA_cor = sum([ crc(count, eff, g_id+'_'+i) for i in ES[event]['CiA'] if g_id+'_'+i in count ])
			ACj_cor = sum([ crc(count, eff, g_id+'_'+i) for i in ES[event]['ACj'] if g_id+'_'+i in count ])
			EXC_cor = sum([ crc(count, eff, g_id+'_'+i) for i in ES[event]['EXC'] if g_id+'_'+i in count ])

			if any([ EXC_raw >= 10, CiA_raw >= 10 and ACj_raw >= 5, CiA_raw >= 5 and ACj_raw >= 10 ]):
				# Order coordinates from ES event
				s, e = natsorted(map(int, event.split('-')))
				try: 
					PSI = 100 * ( (CiA_cor+ACj_cor) / ( (CiA_cor+ACj_cor) + 2 * EXC_cor ) )
					print('{}\tES\t{}-{}\t{:.2f}\t{},{},{},S@{:.2f},{:.2f}\t{}'.format(g_id, s, e, PSI, \
						QC([CiA_raw, ACj_raw], EXC_raw, "ES"), \
						QC([CiA_cor, ACj_cor], EXC_cor, "ES"), \
						RIS([CiA_cor, ACj_cor]), \
						sum([CiA_cor, ACj_cor]), EXC_cor,
						e-s)
					)
				except ZeroDivisionError: pass

def filter_SS_EEJ(eej_iso1, eej_iso2):

	try: 
		rm = []
		for eej1 in eej_iso1:
			s, e = map(int, eej1.split('_')[1].split('-'))

			matches = []
			for eej2 in eej_iso2:
				x, y = map(int, eej2.split('_')[1].split('-'))
				matches.append(is_overlapping(s, e, x, y))

			if not all(matches): rm.append(eej1)

		return [ eej for eej in eej_iso1 if eej not in rm ]
	except ValueError:
		return []

def SS(fasta_files, eff_file, count):
	''' Calculate the "percent splice-site usage" (PSU) 

		Iso1 => ############----------#######
		Iso2 => #######---------------#######

		PSU_iso1 = #Iso1 / (#Iso1 + #Iso2)
		PSU_iso2 = #Iso2 / (#Iso1 + #Iso2)
		#Iso1, #Iso2 => Corrected read counts
	'''
	
	import itertools

	# Parse transcripts
	seqs = {}
	for f in fasta_files:
		seqs = parse_fasta_file(seqs, f)

	# Parse eff file
	eff = parse_eff_count(eff_file)
	g_eff = {}
	for line in open(eff_file):
		eej = line.split('\t')[0]
		g_id = eej.split('_')[0]

		try: g_eff[g_id].append(eej)
		except KeyError: g_eff[g_id] = [eej]

	eprint('#GENE\tSS\tCOORD\tISO1\tISO2')
	print('#GENE\tSS\tCOORD\tPSU\tQUALITY_SCORES\tEVENT_SIZE')

	for g_id in natsorted(seqs):
	#for g_id in ["AT3G56860"]:

		exons = seqs[g_id]
		isFwd = exons.pop('isFwd')

		if len(exons) > 1: # Ignore single transcript genes

			# Throw all exons on one big heap and 
			# ignore first and last exons. FOR NOW...
			alt = { e: [ map(int, e.split('-')) ] for e in list(itertools.chain(*[ natsorted(exons[t].keys())[1:-1] for t in exons ])) }

			# Look for overlaps in the remaining exons
			for exon in natsorted(alt):
				s, e = map(int, exon.split('-'))

				for other_exon in alt:
					if exon != other_exon:
						x, y = map(int, other_exon.split('-'))
						if any([ s <= x <= e and e <= y, x <= s and s <= y <= e ]): alt[exon].append([x, y])

			done = set()

			# If there are any overlaps, collect the EEJs and calculate the PSU
			for iso in natsorted(alt):
				if len(alt[iso]) > 1:

					for i in xrange(1, len(alt[iso])):

						# Only analyse "simple" events where one of the splice sites is shared
						try: 
							shrd = set( alt[iso][0] ).intersection( set( alt[iso][i] ) ).pop()
							iso1, iso2 = set( alt[iso][0] ).symmetric_difference( set( alt[iso][i] ) ) # Get the alternative coordinates

							if frozenset(sorted([shrd, iso1, iso2])) not in done:
								done.add(frozenset(sorted([shrd, iso1, iso2])))

								# Get the EEJ identifiers
								eej_iso1, eej_iso2 = [], []
								index = alt[iso][0].index(shrd)
								for eej in g_eff[g_id]:
									s, e = map(int, eej.split('_')[1].split('-'))
									if not isFwd: s, e = e, s
									if index and e == iso1 or not index and s == iso1: eej_iso1.append(eej)
									if index and e == iso2 or not index and s == iso2: eej_iso2.append(eej)

								# Remove EEJs that are not overlapping
								eej_iso1 = filter_SS_EEJ(eej_iso1, eej_iso2)
								eej_iso2 = filter_SS_EEJ(eej_iso2, eej_iso1)

								iso1_raw = sum([ count[e] if e in count else 0 for e in eej_iso1 ])
								iso2_raw = sum([ count[e] if e in count else 0 for e in eej_iso2 ])

								if index: 
									one = "{}-{}".format(iso1, shrd)
									two = "{}-{}".format(iso2, shrd)
								else:
									one = "{}-{}".format(shrd, iso1)
									two = "{}-{}".format(shrd, iso2)

								eprint('{}\t{}\t{}:{}\t{}\t{}'.format(g_id, 'Alt3' if index else 'Alt5', one, two, iso1_raw, iso2_raw))

								# Continue if any EEJs where found for both events
								#if len(eej_iso1) and len(eej_iso2):
								if iso1_raw > 0 and iso2_raw > 0:

									if iso1_raw+iso2_raw >= 10:

										iso1_cor = sum([ crc(count, eff, e) for e in eej_iso1 ])
										iso2_cor = sum([ crc(count, eff, e) for e in eej_iso2 ])

										PSU1 = 100 * ( iso1_cor / ( iso1_cor + iso2_cor ) )
										PSU2 = 100 * ( iso2_cor / ( iso1_cor + iso2_cor ) )

										print('{}\t{}\t{}:{}\t{:.2f}:{:.2f}\t{},{},NA,S@{:.2f},{:.2f}\t{}'.format(g_id, 'Alt3' if index else 'Alt5', 
											one, two,
											PSU1, PSU2,
											QC([iso1_raw, iso2_raw], "", "SS"),
											QC([iso1_cor, iso2_cor], "", "SS"),
											iso1_cor, iso2_cor,
											abs(alt[iso][0][0] - alt[iso][i][0]) + abs(alt[iso][0][1] - alt[iso][i][1]))
										)

						except KeyError: pass

def ME(eff_file, count):

	# Parse microexon effective file
	eff = parse_eff_count(eff_file)

	# Group per microexon, identify C1 and C2
	# INC => AT1G01010_3913_3945-3954_3996.INC
	# EXC => AT1G01010_3919-3966.EXC

	eprint('#GENE\tCOORD\tINC\tEXC')
	print('#GENE\tEVENT\tCOORD\tPSI\tQUALITY_SCORES\tEVENT_SIZE')

	me = {}
	for line in open(eff_file):
		me_id = line.split('\t')[0]

		if 'INC' in me_id:
			g_id, C1, Aij, C2 = me_id.replace('.INC', '').split('_')
			
			if g_id not in me: me[g_id] = {}
			try: me[g_id][Aij].append([C1, C2])
			except KeyError: me[g_id][Aij] = [[C1, C2]]

	# Get INC (CiAijCj) and EXC (CiCj)
	for g_id in natsorted(me):
		for Aij in me[g_id]:
			INC_raw, EXC_raw = 0, 0
			INC_cor, EXC_cor = 0, 0
			for C1, C2 in me[g_id][Aij]:
				C1AijC2 = '{}_{}_{}_{}.INC'.format(g_id, C1, Aij, C2)
				INC_raw += count[C1AijC2] if C1AijC2 in count else 0
				INC_cor += crc(count, eff, C1AijC2)

				C1C2 = '{}_{}-{}.EXC'.format(g_id, C1, C2)
				EXC_raw += count[C1C2] if C1C2 in count else 0
				EXC_cor += crc(count, eff, C1C2)

			eprint('{}\t{}\t{}\t{}'.format(g_id, Aij, INC_raw, EXC_raw))

			# Calculate PSI => 100 * sum(C1AijC2) / sum(C1AijC2) + sum(C1C2)
			if INC_cor and EXC_cor:
				PSI = 100 * ( INC_cor / ( INC_cor + EXC_raw ) )
				#if 10 <= PSI <= 90 and any([ EXC_raw >= 10, INC_raw >= 10 ]):
				if any([ EXC_raw >= 10, INC_raw >= 10]):
					s, e = map(int, Aij.split('-'))

					print('{}\tME\t{}\t{:.2f}\t{},{},NA,S@{:.2f},{:.2f}\t{}'.format(g_id, Aij, PSI, \
						QC([INC_raw], EXC_raw, "ME"), \
						QC([INC_cor], EXC_cor, "ME"), \
						INC_cor, EXC_cor,
						abs(e-s))
					)

def IR(eff_file, count):
	''' Calculate the "percent intron retention" (PIR)

			  E1I				 IE2
			 ------				------
		########===================########
		   E1	\		I1		  /	  E2
		   		 \               /
		   		  ---------------
		   		   	   E1E2

		PIR = 100 * average(#E1I, #IE2) / (#E1E2 + average(#E1I, #IE2))
		#E1E2, #E1I, #IE2 => Corrected read counts
	'''
	
	from scipy.stats import binom_test
	import numpy as np

	eff = parse_eff_count(eff_file)

	# Get all potential retained introns
	events = {}
	for line in open(eff_file):
		event_id = line.split('\t')[0]
		g_id = event_id.split('_')[0]
		intron = event_id.split('_')[-1]

		try: events[g_id].add(intron)
		except KeyError: events[g_id] = set([intron])

	eprint("#GENE\tIR\tE1E2\tE1I\tIE2\tI1\tIS_BALANCE_OK")
	print('#GENE\tIR\tCOORD\tPIR\tQUALITY_SCORES\tEVENT_SIZE')

	for g_id in natsorted(events):
	#for g_id in ['AT3G56860']:
		for ir in natsorted(events[g_id]):
			E1E2 = '{}_EEJ_{}'.format(g_id, ir)
			E1I = '{}_EIJ_{}'.format(g_id, ir)
			IE2 = '{}_IEJ_{}'.format(g_id, ir)
			I1 = '{}_MID_{}'.format(g_id, ir)

			E1E2_raw = count[E1E2] if E1E2 in count else 0
			E1I_raw = count[E1I] if E1I in count else 0
			IE2_raw = count[IE2] if IE2 in count else 0
			I1_raw = count[I1] if I1 in count else 0

			# Continue if EEJ >= 10 OR EIJ >= 10 and IEJ >= 5 OR EIJ >= 5 and IEJ >= 10
			if any([ E1E2_raw >= 10, (E1I_raw >= 10 and IE2_raw >= 5), (E1I_raw >= 5 and IE2_raw >= 10 ) ]):

				# Get corrected counts
				E1E2_cor = crc(count, eff, E1E2)
				E1I_cor = crc(count, eff, E1I)
				IE2_cor = crc(count, eff, IE2)
				I1_cor = crc(count, eff, I1)

				''' 
					Filter IR events based on (median) read coverage and read balance. 
					The goal of the binomial test (balance) was to exclude events in which there is a 
					high imbalance in read counts among the two exon-intron junctions and the intron 
					body sequence. Such imbalances can arise from neighboring alternative 5' and/or 3'
					splice sites or overlapping genes, confound PIR estimates, and lead to the false 
					detection of IR.

					p-value (binomial{	M = min(#E1I, #IE2, #I), 
										N = min(#E1I, #IE2, #I) + max(#E1I, #IE2, #I),
										P = 1/3.5,
										alternative = lower 
					}) >= 0.05 

					M = number of successes
					N = number of trials
					P = probability of success

					scipy.stats.binom_test(M, N, 1/3.5, alternative="less") 
				'''

				if ( np.median([E1I_cor, IE2_cor, I1_cor]) + E1E2_cor > 10 ) and \
				binom_test(
					min(E1I_cor, IE2_cor, I1_cor),
					min(E1I_cor, IE2_cor, I1_cor) + max(E1I_cor, IE2_cor, I1_cor),
					1/3.5,
					alternative="less"
				) >= 0.05:

					eprint('{}\t{}\t{}\t{}\t{}\t{}\tYES'.format(g_id, ir, E1E2_raw, E1I_raw, IE2_raw, I1_raw))

					try:
						
						# Re-order IR coordinates
						s, e = natsorted(map(int, ir.split('-')))

						PIR = 100 * ( ( ( E1I_cor+IE2_cor ) / 2 ) / ( E1E2_cor + ( ( E1I_cor + IE2_cor ) / 2 ) ) )
						print('{}\tIR\t{}-{}\t{:.2f}\t{},{},{}={}={},S@{:.2f},{:.2f}\t{}'.format(g_id, s, e, PIR, \
							QC([E1I_raw, IE2_raw], E1E2_raw, "IR"), \
							QC([E1I_cor, IE2_cor], E1E2_cor, "IR"), \
							E1I_raw, IE2_raw, E1E2_raw, \
							sum([E1I_cor, IE2_cor]), E1E2_cor,
							e-s)
						)

					except ZeroDivisionError: pass # Ignore cases where E1E2 has no unique mappable positions

				else:
					eprint('{}\t{}\t{}\t{}\t{}\t{}\tNO'.format(g_id, ir, E1E2_raw, E1I_raw, IE2_raw, I1_raw))


def run(fasta_files, eff_file, count_file, AS_type):

	eff = parse_eff_count(eff_file)
	count = parse_eff_count(count_file)

	if AS_type == "ES": # Exon skipping
		ES(fasta_files, eff, count)

	elif AS_type == "ME": # Microexons
		ME(eff_file, count)

	elif AS_type == "IR": # Intron retention
		IR(eff_file, count)

	elif AS_type == "SS": # Alt 3' or Alt 5'
		SS(fasta_files, eff_file, count)

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--fasta', nargs='+', help="Comma separated list of fasta files")
	parser.add_argument('-e', '--eff', help="Effective mappable positions", required=True)
	parser.add_argument('-c', '--count', help="EEJ counts", required=True)
	parser.add_argument('-AS', choices=["ES", "mES", "SS", "IR", "ME"], help="", required=True)
	args = parser.parse_args()

	if args.AS in ["ES", "mES", "SS"] and args.fasta is None:
		parser.error("-AS [{}] requires -f/--fasta".format(args.AS))

	run(args.fasta, args.eff, args.count, args.AS)