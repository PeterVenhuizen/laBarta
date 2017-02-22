#!/usr/bin/env python

"""
Summarize the AStalavista annotation as a counts per event list.
"""

from __future__ import division
import collections
import argparse

def summarize_AStalavista(asta_file, topN=None):
	
	count_AS = collections.Counter()
	for line in open(asta_file):
		c, source, feature, s, e, score, strand, frame, attributes = line.rstrip().split('\t')
		attr = {}
		for a in attributes.split(';'):
			if len(a):
				attr_name, attr_value = a.split(' "')
				attr[attr_name.strip()] = attr_value.replace('\"', '')

		count_AS.update([attr['structure']])

	total = sum(count_AS.values())
	for i, event in enumerate(count_AS.most_common(topN)):
		print '{}\t{}\t{}\t{:.3f}'.format(i+1, event[0], count_AS[event[0]], (count_AS[event[0]]/total)*100)

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-a', '--asta', required=True, help="AStalavista asta output GTF file.")
	parser.add_argument('--topN', default=None, type=int, help="Only return the top N results. Default returns everything.")
	args = parser.parse_args()

	summarize_AStalavista(args.asta, args.topN)