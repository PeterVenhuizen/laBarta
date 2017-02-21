#!/usr/bin/env python

"""
	Beautify GTF file. Remove "Chr" from chromosome, add exon_number, unify transcript ids. 
"""

from natsort import natsorted
import argparse
import re

def beautify(gtf_file):

	# Parse GTF file
	gtf = {}
	chr_sub = re.compile(re.escape('chr'), re.IGNORECASE)
	for line in open(gtf_file):
		c, source, feature, s, e, score, strand, phase, attributes = line.rstrip().split('\t')
		
		# Fix chromosome notation
		c = chr_sub.sub('', c)

		if feature == "exon":
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = a.split(' "')
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			# Fix transcript_id, replace underscore ('_') with dot ('.')
			if "transcript_id" in attr:
				t_id = attr["transcript_id"].replace('_', '.')
				attr["transcript_id"] = t_id
				
				try: 
					gtf[t_id][s+'-'+e] = { 'c': c, 's': s, 'e': e }
				except KeyError:
					gtf[t_id] = { 'info': { 'isFwd': True if strand == '+' else False, 
											'source': source, 
											'feature': feature,
											'attributes': attr
										}
								}
					gtf[t_id][s+'-'+e] = { 'c': c, 's': s, 'e': e }

	# Reformat
	for t_id in natsorted(gtf):

		info = gtf[t_id].pop('info')

		# Sort exons
		ekeys = natsorted(gtf[t_id].keys()) if info['isFwd'] else natsorted(gtf[t_id].keys(), reverse=True)

		for i, key in enumerate(ekeys):
			print '{}\t{}\texon\t{}\t{}\t.\t{}\t.\t{} exon_number "{}";'.format(
				gtf[t_id][key]['c'], 
				info['source'], 
				gtf[t_id][key]['s'],
				gtf[t_id][key]['e'],
				'+' if info['isFwd'] else '-',
				' '.join([ '{} "{}";'.format(k, v) for k, v in info['attributes'].iteritems() ]),
				i+1
			)

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="Input GTF file.")
	args = parser.parse_args()

	beautify(args.gtf)