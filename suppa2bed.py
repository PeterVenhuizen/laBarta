#!/usr/bin/env python

"""
	Extract the alternatively spliced exons defined by SUPPA.

	For each IR event at least two bed tracks will be returned:
		1) The coordinates of just the retained intron
		2) The coordinates of all exons containing the retained intron
"""

import argparse
from fileParser import parse_GTF

def suppa2bed(gtf_file, ioe_file, strict=False):

	# Parse GTF
	gtf = parse_GTF(gtf_file, select_feature="exon", get_introns=False)

	# Parse ioe file
	for line in open(ioe_file):
		seqname, gene_id, event_id, inclusion_transcripts, total_transcripts = line.rstrip().split('\t')

		try: 
			# Parse event_id
			event_type, seqname, coordinates, strand = event_id.split(';')[1].split(':')

			if event_type == 'A3':
				# strict => e1-s2:e1-s3, variable => s2:s3
				pass

			elif event_type == 'A5':
				# strict => e2-s3:e1-s3, variable => s2:s3
				pass

			elif event_type in ['IR', 'RI']:
				# strict => s1:e1-s2:e2, variable => e1:s2
				
				if strict:
					pass
				else: 
					e1, s2 = map(int, coordinates.split('-'))

					# Output the coordinates of just the retained intron
					print '{}\t{}\t{}\t{}:{}\t1000\t{}'.format(seqname, e1, s2, seqname, coordinates, strand)

					# Find exons containing the retained intron
					matches = {}
					for t_id in inclusion_transcripts.split(','):
						for s, e in gtf[t_id]['exons']:
							if s <= e1 <= e and s <= s2 <= e: 
								try: matches['{}-{}'.format(s, e)].append(t_id)
								except KeyError: matches['{}-{}'.format(s, e)] = [t_id]

					for m in matches:
						s, e = m.split('-')
						print '{}\t{}\t{}\t{}:{}:{}\t1000\t{}'.format(seqname, s, e, seqname, coordinates, ','.join(matches[m]), strand)

			elif event_type in ['ES', 'SE']:
				# strict = variable => e1-s2-e2-s3
				pass

			elif event_type == 'MX':
				# strict = variable => e1-s2:e2-s4:e1-s3:e3-s4
				pass

		except IndexError: pass

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--gtf', required=True, help="Genome GTF annotation.")
	parser.add_argument('--ioe', required=True, help="SUPPA events ioe file.")
	parser.add_argument('--strict', choices=[True, False], default=False, help="Were the SUPPA files generated with strict (-b S) or variable (-b V) parameters?")
	args = parser.parse_args()

	suppa2bed(args.gtf, args.ioe, args.strict)