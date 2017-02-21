#!/usr/bin/env python

""" 
Get all constitutive exons from an GTF file. An exon is 
constitutive if it is present in all transcripts of a gene. 
The constitutive exons are output in bed format. 
"""

import argparse
from fileParser import parse_GTF
from natsort import natsorted

def get_constitutive_exons(gtf_file, select_feature="exon", t_id_attr="transcript_id", attr_sep=" "):


	from fileParser import parse_GTF
	gtf = parse_GTF(gtf_file, select_feature, t_id_attr, attr_sep, False)

	# Group per gene
	groups = {}
	for t in gtf:
		g_id = t.split('.')[0]
		try: groups[g_id].append(t)
		except KeyError: groups[g_id] = [t]

	for g_id in natsorted(groups):

		if len(groups[g_id]) > 1:
			trs = groups[g_id]
			first = trs.pop(0)

			for s, e in gtf[first]['exons']:
				if all([ [s, e] in gtf[t_id]['exons'] for t_id in trs ]):
					print '{0}\t{1}\t{2}\t{3}_{0}:{1}-{2}\t1000\t{4}'.format(gtf[t_id]['chr'], s, e, g_id, gtf[t_id]['strand'])

		else:
			t_id = groups[g_id][0]
			for s, e in gtf[t_id]['exons']:
				print '{0}\t{1}\t{2}\t{3}_{0}:{1}-{2}\t1000\t{4}'.format(gtf[t_id]['chr'], s, e, g_id, gtf[t_id]['strand'])

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', default='-', help="Input GTF file")
	parser.add_argument('--select-feature', default="exon", help="Feature to select in the GTF file")
	parser.add_argument('--transcript-attr', default='transcript_id', help="The transcript identifier attribute name")
	parser.add_argument('--attr-sep', default=' "', help="GFF attribute and attribute value separator")
	args = parser.parse_args()

	get_constitutive_exons(args.gtf, args.select_feature, args.transcript_attr, args.attr_sep)