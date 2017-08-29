#!/usr/bin/env python

"""
snpFunctionalAnnotation.py
"""

import argparse
from natsort import natsorted
from sequence import Sequence
from fileParser import parse_GTF

def get_bins(s, e, bin_size=10000):

	b1 = int(int(s)/bin_size)
	b2 = int(int(e)/bin_size)
	bin_code = '{}-{}'.format(b1*bin_size, (b1+1)*bin_size)
	if b1 == b2: return [bin_code]
	else: return [ bin_code, '{}-{}'.format(b2*bin_size, (b2+1)*bin_size) ]

def parse_fasta_file(fasta_file):
	
	from fileParser import yield_fasta

	d = {}
	for seq_record in yield_fasta(fasta_file):
		
		t_id = seq_record.id.split('|')[0]
		locus = seq_record.id.split(' ')[4].split(':')[1]
		start, end = map(int, locus.split('-'))

		try: d[t_id][locus] = seq_record.seq
		except KeyError: d[t_id] = { locus: seq_record.seq }
	
	return d

def snpFunctionalAnnotation(CDS_fasta, CDS_gtf, transcript_gtf, vcf_files):
	
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

	### General annotation (exonic, intronic, etc.) ###
	# Parse transcript GTF
	trs_gtf = parse_GTF(transcript_gtf, get_introns=True)

	typing = { 'E': 'exonic_variant', 'I': 'intronic_variant', '5': 'splice_donor_region_variant', 'F': 'splice_donor_variant', '3': 'splice_acceptor_region_variant', 'T': 'splice_acceptor_variant' }
	intersect = {}
	for t_id in trs_gtf:

		c = trs_gtf[t_id]['chr']
		isFwd = trs_gtf[t_id]['strand'] == '+'
		start = min([ s for s, e in trs_gtf[t_id]['exons'] ])
		end = max([ e for s, e in trs_gtf[t_id]['exons'] ])

		txt_transcript = list('I' * abs(end-start))
		for s, e in trs_gtf[t_id]['exons']: txt_transcript[ s-start:e-start ] = list('E' * (e-s))
		for s, e in trs_gtf[t_id]['introns']:

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
						for i in xrange(pos, pos+(max([len(variants[e]['REF'])] + map(len, variants[e]['ALT'].split(',')) ))):
							try: ann.add(txt_transcript[i-start-1])
							except IndexError: pass

						try: intersect[t_id][e] = { 'label': [ typing[a] for a in ann ], 'effect': '' }
						except KeyError: intersect[t_id] = { e: { 'label': [ typing[a] for a in ann ], 'effect': '' } }

			except KeyError: pass

	### Functional annotation (sense, missense, etc.) ###
	# Parse CDS GTF
	cds_gtf = {}
	for line in open(CDS_gtf):
		if not line.startswith('#'):
			c, source, feature, start, end, score, strand, frame, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = a.split(' "')
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			if feature == "CDS":
				try: cds_gtf[attr["transcript_id"]]["CDS"].append([ int(start), int(end), int(frame) ])
				except KeyError: cds_gtf[attr["transcript_id"]] = { "chr": c, "strand": strand, "CDS": [[ int(start), int(end), int(frame) ]] }

	fa = parse_fasta_file(CDS_fasta)

	#print bins

	#for t_id in natsorted(cds_gtf):
	for t_id in cds_gtf:

		#print t_id
		c = cds_gtf[t_id]['chr']
		isFwd = cds_gtf[t_id]['strand'] == '+' 
		start = min([ s for s, e, p in cds_gtf[t_id]['CDS'] ])
		end = max([ e for s, e, p in cds_gtf[t_id]['CDS'] ])

		start_bin = int(get_bins(start, start)[0].split('-')[0])
		try: end_bin = int(get_bins(end, end)[1].split('-')[1])
		except IndexError: end_bin = int(get_bins(end, end)[0].split('-')[1])

		#print start, end, start_bin, end_bin

		for i in xrange(start_bin, end_bin, 10000):
			
			try:
				for snp in bins[c]['{}-{}'.format(i, i+10000)]:
					chrom, pos, var = snp.split(':')
					ref, alt = var.split('-')
					pos = int(pos)

					# Check if the SNP falls within the CDS coordinates
					if start < pos < end:

						#print t_id, snp, start, end

						# Find the overlapping CDS exon
						found = False
						for j, cds in enumerate(cds_gtf[t_id]['CDS']):
							s, e, phase = cds
							if s < pos < e and not found:

								found = True

								# Only do functional annotation for SNPs, ignore indels
								if len(ref) == 1 and len(alt) == 1:

									# Get the correct codon
									if j: 
										codons = [ fa[t_id]['{}-{}'.format(s, e)][k:k+3] for k in xrange(phase, len(fa[t_id]['{}-{}'.format(s, e)]), 3) ]
										Nth_codon = ((pos-(s+phase)) / 3) if isFwd else (e-(pos+phase)) / 3
										Nth_base = (pos-(s+phase)) % 3 if isFwd else (e-(pos+phase)) % 3
										#print t_id, snp, isFwd, Nth_codon, codons[Nth_codon], Nth_base, codons[Nth_base]
									else: 
										codons = [ fa[t_id]['{}-{}'.format(s, e)][k:k+3] for k in xrange(0, len(fa[t_id]['{}-{}'.format(s, e)]), 3) ]
										Nth_codon = (pos-s) / 3 if isFwd else (e-pos) / 3
										Nth_base = (pos-s) % 3 if isFwd else (e-pos) % 3
										#print t_id, snp, isFwd, Nth_codon, codons[Nth_codon], Nth_base, codons[Nth_codon][Nth_base]

									#print codons, Nth_codon, Nth_base

									# Piece together the codon spanning two exons
									if Nth_codon == -1:
										original = fa[t_id]['{}-{}'.format(last_s, last_e)][-(3-phase):] + fa[t_id]['{}-{}'.format(s, e)][:phase]
										mutated = list(original)

									else:
										original = codons[Nth_codon]
										mutated = list(codons[Nth_codon])

									#if len(original) != 3: print snp

									if isFwd: mutated[Nth_base:Nth_base+1] = alt
									else: mutated[Nth_base:Nth_base+1] = Sequence(alt).get_reverse_complement()
									mutated = ''.join(mutated)

									oAA = Sequence(original).translate()[0]
									mAA = Sequence(mutated).translate()[0]

									# Label the mutation
									if oAA == mAA:
										label = "sense"
									else:
										if oAA == "*": label = "stop-loss"
										elif mAA == "*": label = "stop-gain"
										else: label = "missense"

									#if j and isFwd:
									# Replace exonic_variant with protein_coding, because protein_coding trumps the exonic_variant annotation
									try: intersect[t_id][snp]['label'][intersect[t_id][snp]['label'].index('exonic_variant')] = 'protein_coding'
									except ValueError: intersect[t_id][snp]['label'].append('protein_coding')
									intersect[t_id][snp]['effect'] = '{};{}>{};{}>{}'.format(label, original, mutated, oAA, mAA)
									#print '{}\t{}\tprotein_coding\t{}\t{}>{}\t{}>{}'.format(t_id, snp, label, original, mutated, oAA, mAA)

									#raw_input()

									# Store start and end for lookup if we need this information
									# in case a SNP is in a codon spanning two exons
									last_s, last_e = s, e

								# Indels, just annotate that they are protein coding
								else: 
									intersect[t_id][snp]['label'].append('protein_coding')

			except KeyError: pass

	# Output everything
	for t_id in natsorted(intersect):
		for snp in intersect[t_id]:
			print '{}\t{}\t{}\t{}'.format(t_id, snp, '&'.join(intersect[t_id][snp]['label']), intersect[t_id][snp]['effect'])


if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-g', '--gtf', required=True, help="AtRTD2 Mark's protein CDS GTF file.")
	parser.add_argument('-c', '--cds', required=True, help="AtRTD2 Mark's protein CDS fasta file.")
	parser.add_argument('-t', '--transcripts', required=True, help="AtRTD2 transcriptome GTF file.")
	parser.add_argument('-v', '--vcf', required=True, nargs='+', help="Variant vcf file(s).")
	args = parser.parse_args()

	snpFunctionalAnnotation(args.cds, args.gtf, args.transcripts, args.vcf)