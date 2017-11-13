#!/usr/bin/env python

"""
motifMap.py => Look at the localization of a motif

To-do:
	* Motif mapping to genome
	* GTF feature selection
	* bedtools intersection
	* Genic vs intergenic
	* Coding vs non-coding (intron, exon, UTRs)
	* Proximity of motif to 5' and 3' splice sites
	* GO enrichment
"""

from __future__ import division
from natsort import natsorted
from fileParser import yield_fasta, parse_GTF
from sequence import Sequence
from collections import Counter

import os
import re
import argparse
import subprocess

# http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html
# bedtools intersectBed -a <GTF> -b <MOTIFS.BED> -wo

def map2genome(genome_fasta, motifs, output_dir):
	""" Get all motif hits for each motif and 
	combine the output of all motifs in an extra file. """

	if not os.path.exists(output_dir): os.makedirs(output_dir)

	for motif in motifs:
	
		# Regex for overlapping motif occurences
		p = re.compile('(?=({}))'.format(motif))

		hits = set()
		for record in yield_fasta(genome_fasta):

			# Forward search
			for m in p.finditer(record.seq):
				hits.add((record.id, m.start(), m.start()+len(motif), '+'))

			# Reverse search
			seq = Sequence(record.seq).get_reverse_complement()
			seq_len = len(record.seq)
			for m in p.finditer(seq):
				hits.add((record.id, seq_len-(m.start()+len(motif)), seq_len-m.start(), '-'))

		with open('{}/{}.bed'.format(output_dir, motif), 'w') as fout: 
			for hit in natsorted(hits):
				c, s, e, strand = hit
				fout.write( '{}\t{}\t{}\t{}\t1000\t{}\n'.format(c, s, e, motif, strand) )

	# Get all motif hits in one file
	subprocess.call("cat {0}/*bed > {0}/all_motif_hits.bed".format(output_dir), shell=True)

def GxIG(motif_hit_files, gtf_file, output_dir):
	""" Summarize the genic and intergenic hits based on 
	the provided GTF file. Also split the genic hits in 
	strand-specific and non-specific hits. """

	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Create the "gene" bed file, so that the intronic hits
	# are also included as genic hits.
	gtf = parse_GTF(gtf_file)
	genes = {}
	for t_id in gtf:
		g_id = t_id.split('.')[0]
		start = min([ x[0] for x in gtf[t_id]['exons']])
		end = max([ x[1] for x in gtf[t_id]['exons']])
		try:
			genes[g_id]['start'] = start if start < genes[g_id]['start'] else genes[g_id]['start']
			genes[g_id]['end'] = end if end > genes[g_id]['end'] else genes[g_id]['end']
		except KeyError: genes[g_id] = { 'chr': gtf[t_id]['chr'], 'strand': gtf[t_id]['strand'], 'start': start, 'end': end }

	with open("{}/genes.bed".format(output_dir), 'w') as fout:
		for g_id in natsorted(genes):
			fout.write( '{}\t{}\t{}\t{}\t1000\t{}\n'.format(genes[g_id]['chr'], genes[g_id]['start'], genes[g_id]['end'], g_id, genes[g_id]['strand']) )

	with open("{}/GxIG.txt".format(output_dir), 'w') as fout: 
		fout.write('MOTIF\tTOTAL_HITS\tGENIC\tINTERGENIC\tON_STRAND\tOPPOSITE_STRAND\n')
		for motif_file in motif_hit_files:

			motif = motif_file.split('/')[-1].replace('.bed', '')

			# Get the number of motif hits
			n_hits = sum([ 1 for line in open(motif_file) ])

			# Get all genic hits
			p = subprocess.Popen(['intersectBed -a {} -b {}/genes.bed -u | wc -l'.format(motif_file, output_dir)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			genic_hits = int(out.rstrip())

			# Get all intergenic hits
			p = subprocess.Popen(['intersectBed -a {} -b {}/genes.bed -v | wc -l'.format(motif_file, output_dir)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			intergenic_hits = int(out.rstrip())

			# Get the genes with non strand specific hits
			subprocess.call('intersectBed -a {0}/genes.bed -b {1} -u | awk \'{{ print $4 }}\' | sort | uniq > {0}/{2}_genic_non_specific.txt'.format(output_dir, motif_file, motif), shell=True)

			# Get strand-specific genic hits
			p = subprocess.Popen(['intersectBed -a {} -b {}/genes.bed -s -u | wc -l'.format(motif_file, output_dir)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			ss_hits = int(out.rstrip())

			# Get opposite strand genic hits
			p = subprocess.Popen(['intersectBed -a {} -b {}/genes.bed -S -u | wc -l'.format(motif_file, output_dir)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			fs_hits = int(out.rstrip())

			# Get the genes with strand specific hits
			subprocess.call('intersectBed -a {0}/genes.bed -b {1} -s -u | awk \'{{ print $4 }}\' | sort | uniq > {0}/{2}_genic_specific.txt'.format(output_dir, motif_file, motif), shell=True)

			fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(motif, n_hits, genic_hits, intergenic_hits, ss_hits, fs_hits))

def ExI(gtf_file, motif_hit_files, reference_list, output_dir):
	"""
		Get the motif counts for the exons and introns
		of the reference transcripts. 
	"""

	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Create the intron/exon bed files for 
	# the reference isoforms

	ref = set([ line.rstrip() for line in open(reference_list) ])
	exonOut = open('{}/ref_exons.bed'.format(output_dir), 'w')
	intronOut = open('{}/ref_introns.bed'.format(output_dir), 'w')

	gtf = parse_GTF(gtf_file)
	for t_id in natsorted(gtf):
		if t_id in ref:

			c, strand, exons = [gtf[t_id][x] for x in ['chr', 'strand', 'exons']]

			# Output exons
			for i, x in enumerate(exons): 
				exonOut.write( '{}\t{}\t{}\t{}:{}\t1000\t{}\n'.format(c, x[0]-1, x[1], t_id, i+1, strand) )

			# Output introns
			for i in xrange(len(exons)-1):
				x1 = exons[i]
				x2 = exons[i+1]
				intron = [x1[1]+1, x2[0]-1] if strand == '+' else [x2[1]+1, x1[0]-1]
				intronOut.write( '{}\t{}\t{}\t{}:{}\t1000\t{}\n'.format(c, intron[0]-1, intron[1], t_id, i+1, strand) )

	exonOut.close(), intronOut.close()

	with open("{}/ExI.txt".format(output_dir), 'w') as fout:
		fout.write( 'MOTIF\tTOTAL_HITS\tEXONIC\tINTRONIC\n' )
		for motif_file in motif_hit_files:

			motif = motif_file.split('/')[-1].replace('.bed', '')

			# Exonic on-strand hits
			p = subprocess.Popen(['intersectBed -a {} -b {}/ref_exons.bed -f 1 -u -s | wc -l'.format(motif_file, output_dir)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			exonic_hits = int(out.rstrip())

			# Intronic on-strand hits
			p = subprocess.Popen(['intersectBed -a {} -b {}/ref_introns.bed -f 1 -u -s | wc -l'.format(motif_file, output_dir)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			intronic_hits = int(out.rstrip())

			fout.write( '{}\t{}\t{}\t{}\n'.format(motif, exonic_hits+intronic_hits, exonic_hits, intronic_hits) )

def motifDensity(gtf_file, motif_hit_files, reference_list, output_dir, create_bins=False):
	""" 
		Get the motif counts for specific bins. The
		bins are as follows:
			* Upstream exon midpoint to 30bp upstream of exon end
			* 30bp upstream of exon end to 5' splice site
			* 14bp around 5' splice site
			* 30bp downstream of the 5' splice site
			* Everything between 30bp downstream and 30bp upstream of the intron
			* 30bp upstream of the 3' splice site
			* 14bp around 3' splice site
			* 30bp downstream of exon start
			* 30bp from start to midpoint
	"""

	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Create the bins
	if create_bins:
		gtf = parse_GTF(gtf_file)

		ref = set([ line.rstrip() for line in open(reference_list) ])

		with open('{}/gtf2bins.bed'.format(output_dir), 'w') as fout: 
			for t_id in gtf:
				if t_id in ref:
					c = gtf[t_id]['chr']
					exons = gtf[t_id]['exons']
					strand = gtf[t_id]['strand']
					if len(exons) >= 2:

						for i in xrange(len(exons)-1):

							ex1 = exons[i]
							ex2 = exons[i+1]
							intron = [ex1[1]+1, ex2[0]-1] if strand == '+' else [ex2[1]+1, ex1[0]-1]
							if strand == '-': ex1, ex2 = ex2, ex1

							# 5' splice site
							ex1_len = abs(ex1[1]-ex1[0])
							ex1_mid = int(ex1[0] + round(ex1_len / 2))
							if ex1_len > 74:
								fout.write( '{}\t{}\t{}\tmid_30bp\t1000\t{}\n'.format(c, ex1_mid, ex1[1]-31, strand) )
								fout.write( '{}\t{}\t{}\t30bp_{}prime\t1000\t{}\n'.format(c, ex1[1]-30, ex1[1], '5' if strand == '+' else '3', strand) )
							elif ex1_len > 15:
								fout.write( '{}\t{}\t{}\t30bp_{}prime\t1000\t{}\n'.format(c, ex1_mid, ex1[1], '5' if strand == '+' else '3', strand) )
							fout.write( '{}\t{}\t{}\t{}prime\t1000\t{}\n'.format(c, ex1[1]-6, ex1[1]+8, '5' if strand == '+' else '3', strand) )

							# intron
							intron_len = abs(intron[1]-intron[0])
							if intron_len > 74:
								fout.write( '{}\t{}\t{}\t{}prime_30bp\t1000\t{}\n'.format(c, intron[0]+1, intron[0]+30, '5' if strand == '+' else '3', strand) )
								x = natsorted([intron[0], intron[1]])
								fout.write( '{}\t{}\t{}\tintron_mid\t1000\t{}\n'.format(c, x[0]+30, x[1]-30, strand) )
								fout.write( '{}\t{}\t{}\t30bp_{}prime\t1000\t{}\n'.format(c, intron[1]-31, intron[1]-2, '3' if strand == '+' else '5', strand) )
							else:
								fout.write( '{}\t{}\t{}\t{}prime_30bp\t1000\t{}\n'.format(c, intron[0]+1, int(intron[0]+round(intron_len / 2)), '5' if strand == '+' else '3', strand) )
								fout.write( '{}\t{}\t{}\t30bp_{}prime\t1000\t{}\n'.format(c, int(intron[1]-round(intron_len / 2)), intron[1]-2, '3' if strand == '+' else '5', strand) )

							# 3' splice site
							ex2_len = abs(ex2[1]-ex2[0])
							ex2_mid = int(ex2[0] + round(ex2_len / 2))
							fout.write( '{}\t{}\t{}\t{}prime\t1000\t{}\n'.format(c, ex2[0]-9, ex2[0]+5, '3' if strand == '+' else '5', strand) )
							if ex2_len > 74:
								fout.write( '{}\t{}\t{}\t{}prime_30bp\t1000\t{}\n'.format(c, ex2[0], ex2[0]+30, '3' if strand == '+' else '5', strand) )
								fout.write( '{}\t{}\t{}\t30bp_mid\t1000\t{}\n'.format(c, ex2[0]+31, ex2_mid, strand) )
							elif ex2_len > 15:
								fout.write( '{}\t{}\t{}\t{}prime_30bp\t1000\t{}\n'.format(c, ex2[0], ex2_mid, '3' if strand == '+' else '5', strand) )

	# Get the total length of each bin
	scale = Counter()
	for line in open('{}/gtf2bins.bed'.format(output_dir)):
		c, s, e, b = line.split('\t')[:4]
		scale[b] += abs(int(e)-int(s))

	# Get the motif hits
	with open('{}/motifDensity.txt'.format(output_dir), 'w') as fout:
		fout.write('MOTIF\tSTRANDNESS\tMID_30BP\t30BP_5PRIME\t5PRIME\t5PRIME_30BP\tINTRON_MID\t30BP_3PRIME\t3PRIME\t3PRIME_30BP\t30BP_MID\n')
		for motif_file in motif_hit_files:
			motif = motif_file.split('/')[-1].replace('.bed', '')

			# on-strand
			on_strand = Counter()
			subprocess.call("intersectBed -a {0}/gtf2bins.bed -b {1} -s -F 1 -wo > {0}/on-strand.txt".format(output_dir, motif_file), shell=True)
			for line in open('{}/on-strand.txt'.format(output_dir)):
				on_strand[line.split('\t')[3]] += 1
			fout.write( '{}\ton-strand\t{}\t{}\n'.format(motif, sum(on_strand.values()), '\t'.join([ '{:.6f}'.format(on_strand[x] * (1000/scale[x])) for x in ['mid_30bp', '30bp_5prime', '5prime', '5prime_30bp', 'intron_mid', '30bp_3prime', '3prime', '3prime_30bp', '30bp_mid'] ])) )

			# off-strand
			off_strand = Counter()
			subprocess.call("intersectBed -a {0}/gtf2bins.bed -b {1} -S -F 1 -wo > {0}/off-strand.txt".format(output_dir, motif_file), shell=True)
			for line in open('{}/off-strand.txt'.format(output_dir)):
				off_strand[line.split('\t')[3]] += 1
			fout.write( '{}\toff-strand\t{}\t{}\n'.format(motif, sum(off_strand.values()), '\t'.join([ '{:.6f}'.format(off_strand[x] * (1000/scale[x])) for x in ['mid_30bp', '30bp_5prime', '5prime', '5prime_30bp', 'intron_mid', '30bp_3prime', '3prime', '3prime_30bp', '30bp_mid'] ])) )

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-d', '--dir', required=True, help="Output directory.")
	subparsers = parser.add_subparsers(dest='command', help="Sub-command help.")

	parser_a = subparsers.add_parser('map2genome', help="Map motif(s) to the genomic fasta sequence.")
	parser_a.add_argument('--genome', required=True, help="Genome fasta sequence.")
	parser_a.add_argument('--motif', nargs='+', required=True, help="List of motifs.")

	parser_b = subparsers.add_parser('GxIG', help="Count genic vs intergenic hits.")
	parser_b.add_argument('--motif-hits', nargs='+', required=True, help="List of motif hit files (from map2genome).")
	parser_b.add_argument('-g', '--gtf', required=True, help="Transcript GTF annotation.")

	parser_c = subparsers.add_parser('ExI', help="Count the exonic vs intronic hits.")
	parser_c.add_argument('--motif-hits', nargs='+', required=True, help="List of motif hit files (from map2genome).")
	parser_c.add_argument('-g', '--gtf', required=True, help="Transcript GTF annotation.")
	parser_c.add_argument('-r', '--ref', required=True, help="Representative gene models.")
	
	parser_d = subparsers.add_parser('motifDensity', help="Get the motif density.")
	parser_d.add_argument('--motif-hits', nargs='+', required=True, help="List of motif hit files (from map2genome).")
	parser_d.add_argument('-g', '--gtf', required=True, help="Transcript GTF annotation.")
	parser_d.add_argument('-r', '--ref', required=True, help="Representative gene models.")
	parser_d.add_argument('--create-bins', action='store_true')

	args = parser.parse_args()

	if args.command == 'map2genome':
		map2genome(args.genome, args.motif, args.dir)
	elif args.command == "GxIG":
		GxIG(args.motif_hits, args.gtf, args.dir)
	elif args.command == "ExI":
		ExI(args.gtf, args.motif_hits, args.ref, args.dir)
	elif args.command == "motifDensity":
		motifDensity(args.gtf, args.motif_hits, args.ref, args.dir, args.create_bins)