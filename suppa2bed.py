#!/usr/bin/env python

"""
	Extract the alternatively spliced exons defined by SUPPA.

	For each IR event at least two bed tracks will be returned:
		1) The coordinates of just the retained intron
		2) The coordinates of all exons containing the retained intron
"""

from __future__ import print_function
import argparse
import sys
from fileParser import parse_GTF

def eprint(*args, **kwargs):
	print(*args, file=sys.stderr, **kwargs)

def suppa2bed(gtf_file, ioe_file, strict="F"):

	strict = True if strict == "T" else False

	# Parse GTF
	gtf = parse_GTF(gtf_file, select_feature="exon", get_introns=False)

	# Parse ioe file
	for line in open(ioe_file):
		seqname, gene_id, event_id, inclusion_transcripts, total_transcripts = line.rstrip().split('\t')

		try: 
			# Parse event_id
			event_type = event_id.split(';')[1].split(':')[0]

			if event_type == 'A3':
				# strict => e1-s2:e1-s3, variable => s2:s3
				
				if strict:
					event_type, seqname, e1_s2, e1_s3, strand = event_id.split(';')[1].split(':')
					e1, s2 = e1_s2.split('-')
					print('{0}\t{1}\t{2}\t{0}:{1}-{2}\t1000\t{3}'.format(seqname, e1, s2, strand))

					e1, s3 = e1_s3.split('-')
					print('{0}\t{1}\t{2}\t{0}:{1}-{2}\t1000\t{3}'.format(seqname, e1, s3, strand))

				else: 
					pass

			elif event_type == 'A5':
				# strict => e2-s3:e1-s3, variable => s2:s3
				if strict:
					event_type, seqname, e2_s3, e1_s3, strand = event_id.split(';')[1].split(':')
					e2, s3 = e2_s3.split('-')
					print('{0}\t{1}\t{2}\t{0}:{1}-{2}\t1000\t{3}'.format(seqname, e2, s3, strand))

					e1, s3 = e1_s3.split('-')
					print('{0}\t{1}\t{2}\t{0}:{1}-{2}\t1000\t{3}'.format(seqname, e1, s3, strand))

				else: 
					pass

			elif event_type in ['IR', 'RI']:
				# strict => s1:e1-s2:e2, variable => e1:s2
				
				if strict:
					event_type, seqname, s1, e1_s2, e2, strand = event_id.split(';')[1].split(':')
					e1, s2 = e1_s2.split('-')
					coordinates = '{}:{}:{}'.format(s1, e1_s2, e2)

					# Output the coordinates of just the retained intron
					print('{}\t{}\t{}\t{}:{}\t1000\t{}'.format(seqname, e1, s2, seqname, coordinates, strand))

					# Output IR exon
					eprint('{}\t{}\t{}\t{}:{}:{}\t1000\t{}'.format(seqname, s1, e2, seqname, coordinates, inclusion_transcripts, strand))
				else: 

					event_type, seqname, coordinates, strand = event_id.split(';')[1].split(':')
					e1, s2 = map(int, coordinates.split('-'))

					# Output the coordinates of just the retained intron
					print('{}\t{}\t{}\t{}:{}\t1000\t{}'.format(seqname, e1, s2, seqname, coordinates, strand))

					# Find exons containing the retained intron
					matches = {}
					for t_id in inclusion_transcripts.split(','):
						for s, e in gtf[t_id]['exons']:
							if s <= e1 <= e and s <= s2 <= e: 
								try: matches['{}-{}'.format(s, e)].append(t_id)
								except KeyError: matches['{}-{}'.format(s, e)] = [t_id]

					for m in matches:
						s, e = m.split('-')
						eprint('{}\t{}\t{}\t{}:{}:{}\t1000\t{}'.format(seqname, s, e, seqname, coordinates, ','.join(matches[m]), strand))

			elif event_type in ['ES', 'SE']:
				# strict = variable => e1-s2-e2-s3
								
				if strict:
					event_type, seqname, e1_s2, e2_s3, strand = event_id.split(';')[1].split(':')
					e1, s2 = e1_s2.split('-')
					e2, s3 = e2_s3.split('-')

					print('{0}\t{1}\t{2}\t{0}:{1}-{2}\t1000\t{3}'.format(seqname, e1, s3, strand))

				else: 
					pass

			elif event_type == 'MX':
				# strict = variable => e1-s2:e2-s4:e1-s3:e3-s4
				pass

		except IndexError: pass

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--gtf', required=True, help="Genome GTF annotation.")
	parser.add_argument('--ioe', required=True, help="SUPPA events ioe file.")
	parser.add_argument('--strict', choices=["T", "F"], default="F", help="Were the SUPPA files generated with strict (-b S) or variable (-b V) parameters?")
	args = parser.parse_args()

	suppa2bed(args.gtf, args.ioe, args.strict)