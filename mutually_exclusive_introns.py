#!/usr/bin/env python

"""
mxi.py => Detect Mutually Exclusive Introns in a given GTF file. 

	Okay scenario 1
	---|||||-----|||||||||||||||---
	---|||||||||||||||-----|||||---
	---|||||-----|||||-----|||||---
	
	Okay scenario 2
	---|||||-----|||||-----|||||||||||||||---
	---|||||||||||||||-----|||||-----|||||---
	---|||||-----|||||-----|||||-----|||||---

	Bad scenario 1
	---|||||||||||||||||||||||||---
	---|||||-----|||||||||||||||---
	---|||||||||||||||-----|||||---
	---|||||-----|||||-----|||||---

	AStalavista MXI code => 1^2-,3^4-

	To-do:
		* Exclude the cases in which all introns are retained
		* Consider distant introns for MXI

"""

import argparse
from fileParser import parse_GTF
from natsort import natsorted

def is_contained_in(s, e, x, y):
	''' Check if s-e is contained in x-y '''
	s, e = natsorted((s, e))
	x, y = natsorted((x, y))
	return x < s < y and x < e < y

def mxi(gtf_file):

	gtf = parse_GTF(gtf_file, select_feature="exon", get_introns=True)

	# Groups per gene
	groups = {}
	for t_id in gtf:
		g_id = t_id.split('.')[0]
		try: groups[g_id].append(t_id)
		except KeyError: groups[g_id] = [t_id]

	print "seqname\tgene_id\tevent_id\talternative_transcripts\ttotal_transcripts"
	for g_id in natsorted(groups):

		trs = set([ t_id for t_id in groups[g_id] ])
		c, strand = [ gtf[next(iter(trs))][x] for x in ['chr', 'strand'] ]

		# Get the retained introns
		# Get the inclusion (containing the intron) transcripts
		introns = {}
		for t_id in trs: 
			for s, e in gtf[t_id]['introns']:
				try: introns['{}-{}'.format(s, e)]['inc'].add(t_id)
				except KeyError: introns['{}-{}'.format(s, e)] = { 'inc': set([t_id]), 'exc': set([]) }

		# Get the exclusion transcripts
		for i in introns:
			s, e = map(int, i.split('-'))
			for t_id in trs.difference(introns[i]['inc']):
				for x, y in gtf[t_id]['exons']:
					if is_contained_in(s, e, x, y):
						introns[i]['exc'].add(t_id)
						break

		mxi = {}
		print gtf['AT2G46610.CR1']['introns']
		print gtf['AT2G46610.CR3']['introns']

		# I don't want this one => AT2G46610;MXI:2:19137681-19137830:19138514-19138640:- 
		# Find a way to avoid this and other cases like it

		for i1 in natsorted(introns):
			s, e = map(int, i1.split('-'))
			#print s, e
			for t1 in introns[i1]['inc']:
				for t2 in introns[i1]['exc']:

					# Get the next intron
					if len(gtf[t2]['introns']) > 0:

						if strand == '+': n = sum([ y < e for x, y in gtf[t2]['introns'] ])
						elif strand == '-': 
							n = sum([ y > s for x, y in gtf[t2]['introns'] ])
							print t2, [ y > s for x, y in gtf[t2]['introns'] ]
							print gtf[t2]['introns']
							n = n-1 if n > 0 else n

						try: 
							#n = 0 if n == 1 and len(gtf[t2]['introns']) == 1 else n

							#print [ y < e for x, y in gtf[t2]['introns'] ], gtf[t2]['introns']
							#print t2, n, gtf[t2]['introns']

							i2 = '{}-{}'.format(*gtf[t2]['introns'][n])
							#print i1, i2
							#raw_input()

							if t1 in introns[i2]['exc'] and t2 in introns[i2]['inc']:
								try:
									mxi['{}:{}'.format(i1, i2)]['i1'].add(t1)
									mxi['{}:{}'.format(i1, i2)]['i2'].add(t2)
								except KeyError:
									if '{}:{}'.format(*natsorted([i1, i2])) not in mxi:
										mxi['{}:{}'.format(i1, i2)] = { 'i1': set([t1]), 'i2': set([t2]) }
						except IndexError: pass

		for x in natsorted(mxi):
			print '{0}\t{1}\t{1};MXI:{0}:{2}:{3}\t{4}\t{5}'.format(c, g_id, x, strand, ','.join(mxi[x]['i1']), ','.join(mxi[x]['i1'] | mxi[x]['i2']))
			#print '{}\t{}\t{}\t{}\t1000\t{}'.format(c, x.split('-')[0], x.split('-')[-1], x, strand)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Transcripts in GTF format.")
	args = parser.parse_args()

	mxi(args.gtf)