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

	with open("{}/GxIG.txt".format(output_dir), 'w') as fout: 
		fout.write('MOTIF\tTOTAL_HITS\tGENIC\tINTERGENIC\tON_STRAND\tOPPOSITE_STRAND\n')
		for motif_file in motif_hit_files:

			motif = motif_file.split('/')[-1].replace('.bed', '')

			# Get the number of motif hits
			n_hits = sum([ 1 for line in open(motif_file) ])

			# Get all genic hits
			p = subprocess.Popen(['intersectBed -a {} -b {} -u | wc -l'.format(motif_file, gtf_file)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			genic_hits = int(out.rstrip())

			# Get the genes with non strand specific hits
			subprocess.call('intersectBed -a {} -b {} -u | python laBarta/get_gene_id.py | sort | uniq > {}/{}_genic_non_specific.txt'.format(gtf_file, motif_file, output_dir, motif), shell=True)

			# Get strand-specific genic hits
			p = subprocess.Popen(['intersectBed -a {} -b {} -s -u | wc -l'.format(motif_file, gtf_file)], stdout=subprocess.PIPE, shell=True)
			out, err = p.communicate()
			ss_hits = int(out.rstrip())
			fs_hits = genic_hits - ss_hits

			# Get the genes with strand specific hits
			subprocess.call('intersectBed -a {} -b {} -s -u | python laBarta/get_gene_id.py | sort | uniq > {}/{}_genic_specific.txt'.format(gtf_file, motif_file, output_dir, motif), shell=True)

			# Get the number of intergenic hits
			intergenic_hits = n_hits - genic_hits

			fout.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(motif, n_hits, genic_hits, intergenic_hits, ss_hits, fs_hits))

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
							intron = [ex1[1]+1, ex2[0]-1] if strand == '+' else [ex1[0]-1, ex2[1]+1]

							# 5' splice site
							ex1_len = abs(ex1[1]-ex1[0])
							ex1_mid = int(ex1[0] + round(ex1_len / 2))
							if ex1_len > 74:
								fout.write( '{}\t{}\t{}\tmid_30bp\t1000\t{}\n'.format(c, ex1_mid, ex1[1]-31, strand) )
								fout.write( '{}\t{}\t{}\t30bp_5prime\t1000\t{}\n'.format(c, ex1[1]-30, ex1[1], strand) )
							elif ex1_len > 15:
								fout.write( '{}\t{}\t{}\t30bp_5prime\t1000\t{}\n'.format(c, ex1_mid, ex1[1], strand) )
							fout.write( '{}\t{}\t{}\t5prime\t1000\t{}\n'.format(c, ex1[1]-5, ex1[1]+8, strand) )

							# intron
							intron_len = abs(intron[1]-intron[0])
							if intron_len > 74:
								fout.write( '{}\t{}\t{}\t5prime_30bp\t1000\t{}\n'.format(c, intron[0]+2, intron[0]+30, strand) )
								x = natsorted([intron[0], intron[1]])
								fout.write( '{}\t{}\t{}\tintron_mid\t1000\t{}\n'.format(c, x[0]+30, x[1]-30, strand) )
								fout.write( '{}\t{}\t{}\t30bp_3prime\t1000\t{}\n'.format(c, intron[1]-30, intron[1]-2, strand) )
							else:
								fout.write( '{}\t{}\t{}\t5prime_30bp\t1000\t{}\n'.format(c, intron[0]+2, int(intron[0]+round(intron_len / 2)), strand) )
								fout.write( '{}\t{}\t{}\t30bp_3prime\t1000\t{}\n'.format(c, int(intron[1]-round(intron_len / 2)), intron[1]-2, strand) )

							# 3' splice site
							ex2_len = abs(ex2[1]-ex2[0])
							ex2_mid = int(ex2[0] + round(ex2_len / 2))
							fout.write( '{}\t{}\t{}\t3prime\t1000\t{}\n'.format(c, ex2[0]-8, ex2[0]+5, strand) )
							if ex2_len > 74:
								fout.write( '{}\t{}\t{}\t3prime_30bp\t1000\t{}\n'.format(c, ex2[0], ex2[0]+30, strand) )
								fout.write( '{}\t{}\t{}\t30bp_mid\t1000\t{}\n'.format(c, ex2[0]+31, ex2_mid, strand) )
							elif ex2_len > 15:
								fout.write( '{}\t{}\t{}\t3prime_30bp\t1000\t{}\n'.format(c, ex2[0], ex2_mid, strand) )

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
	
	parser_c = subparsers.add_parser('motifDensity', help="Get the motif density.")
	parser_c.add_argument('--motif-hits', nargs='+', required=True, help="List of motif hit files (from map2genome).")
	parser_c.add_argument('-g', '--gtf', required=True, help="Transcript GTF annotation.")
	parser_c.add_argument('-r', '--ref', required=True, help="Representative gene models.")
	parser_c.add_argument('--create-bins', action='store_true')

	args = parser.parse_args()

	if args.command == 'map2genome':
		map2genome(args.genome, args.motif, args.dir)
	elif args.command == "GxIG":
		GxIG(args.motif_hits, args.gtf, args.dir)
	elif args.command == "motifDensity":
		motifDensity(args.gtf, args.motif_hits, args.ref, args.dir, args.create_bins)