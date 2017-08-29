#!/usr/bin/env python

"""
Parse GATK/snpEff output and intersect with GTF
"""

import argparse
from natsort import natsorted
from fileParser import parse_GTF

def get_bins(s, e, bin_size=10000):

	b1 = int(int(s)/bin_size)
	b2 = int(int(e)/bin_size)
	bin_code = '{}-{}'.format(b1*bin_size, (b1+1)*bin_size)
	if b1 == b2: return [bin_code]
	else: return [ bin_code, '{}-{}'.format(b2*bin_size, (b2+1)*bin_size) ]

def intersect_vcf(vcf_files, gtf_file):

	# Parse vcf
	variants = {}
	for f in vcf_files:
		for line in open(f):
			if not line.startswith('#'):
				CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO = line.rstrip().split('\t')[:8]
				variants['{}:{}:{}-{}'.format(CHROM, POS, REF, ALT)] = { 'ID': ID, 'REF': REF, 'ALT': ALT, 'QUAL': float(QUAL), 'FILTER': FILTER, 'INFO': INFO }

	# Bin the vcf events
	bins = {}
	for e in natsorted(variants):
		if variants[e]['FILTER'] == 'PASS':
			c, pos, var = e.split(':')

			if c not in bins: bins[c] = {}
			bin_codes = get_bins(int(pos), int(pos))
			for b in bin_codes:
				try: bins[c][b].append(e)
				except KeyError: bins[c][b] = [e]

	# Parse GTF
	gtf = parse_GTF(gtf_file, get_introns=True)

	# Intersect GTF and vcf
	typing = { 'E': 'exonic_variant', 'I': 'intronic_variant', '5': 'splice_donor_region_variant', 'F': 'splice_donor_variant', '3': 'splice_acceptor_region_variant', 'T': 'splice_acceptor_variant' }
	intersect = {}
	for t_id in gtf:

		c = gtf[t_id]['chr']
		isFwd = gtf[t_id]['strand'] == '+'
		start = min([ s for s, e in gtf[t_id]['exons'] ])
		end = max([ e for s, e in gtf[t_id]['exons'] ])

		txt_transcript = list('I' * abs(end-start))
		for s, e in gtf[t_id]['exons']: txt_transcript[ s-start:e-start ] = list('E' * (e-s))
		for s, e in gtf[t_id]['introns']:

			if isFwd:
				txt_transcript[ s-start-3:s-start+10 ] = list('5' * 13)
				txt_transcript[ e-start-14:e-start+3 ] = list('3' * 17)
				txt_transcript[ s-start:s-start+2 ] = list('FF')
				txt_transcript[ e-start-2:e-start ] = list('TT')

			else: 
				txt_transcript[ s-start-3:s-start+14 ] = list('3' * 17)
				txt_transcript[ e-start-10:e-start+3 ] = list('5' * 13)
				txt_transcript[ s-start:s-start+2 ] = list('TT')
				txt_transcript[ e-start-2:e-start ] = list('FF')

		start_bin = int(get_bins(start, start)[0].split('-')[0])
		try: end_bin = int(get_bins(end, end)[1].split('-')[1])
		except IndexError: end_bin = int(get_bins(end, end)[0].split('-')[1])

		for i in xrange(start_bin, end_bin, 10000):
			try: 
				for e in bins[c]['{}-{}'.format(i, i+10000)]:
					chrom, pos, var = e.split(':')
					pos = int(pos)
					if start < pos < end:

						# Variant typing
						ann = set()
						for i in xrange(pos, pos+(max([len(variants[e]['REF'])] + map(len, variants[e]['ALT'].split(','))))):
							try: ann.add(txt_transcript[i-start-1])
							except IndexError: pass
						var_type = '&'.join([ typing[a] for a in ann ])

						try: intersect[t_id].append('{}:{}'.format(e, var_type))
						except KeyError: intersect[t_id] = [ '{}:{}'.format(e, var_type) ]

			except KeyError: pass

	for t_id in natsorted(intersect):
		print '{}\t{}\t{}'.format(t_id, len(intersect[t_id]), ', '.join(intersect[t_id]))


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-v', '--vcf', nargs='+', required=True, help="SNP/Indel vcf file.")
	parser.add_argument('-g', '--gtf', required=True, help="Transcript GTF annotation.")
	args = parser.parse_args()

	intersect_vcf(args.vcf, args.gtf)