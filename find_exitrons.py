#!/usr/bin/env python

'''
	Multi-purpose script for exitron detection. Some functions require samtools
	(tested with v1.2) and/or bedtools (tested with v2.26.0) to be installed. 

	Created on:		16-08-2016 (Peter Venhuizen)
	Last edited: 	22-02-2017 (Peter Venhuizen)
'''

from __future__ import division
from natsort import natsorted
import argparse
import subprocess
import os
import sys
import re

from fileParser import yield_junctions, yield_fasta
from sequence import Sequence

def query_yes_no(question, default="y"):
    """Ask a yes/no question via raw_input() and return their answer.
    "question" is a string that is presented to the user.
    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "no": False, "n": False}

    while True:
        sys.stderr.write(question + "[y/n]")
        choice = raw_input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stderr.write("Please respond with 'yes' or 'no' (or 'y' or 'n').\n")

def add_slash(path):
	return path if path.endswith('/') else path + '/'

def yield_bed(bed_file):

	for line in open(bed_file):
		cols = line.rstrip().split('\t')
		d = { 
			'c': cols[0],
			's': int(cols[1]),
			'e': int(cols[2]),
			'id': cols[3],
			'score': int(cols[4]),
			'strand': cols[5]
		}
		if len(cols) > 6: d['phase'] = int(cols[6])
		yield(d)

def get_shared_junctions(junction_files):
	return set.intersection(*[ set([ (j['chr'], j['junc_start'], j['junc_end'], j['strand']) for j in yield_junctions(f) if j['depth'] >= 3 ]) for f in junction_files ])

def get_bin(s, e, bin_size=5000):

	bin1 = int(int(s)/bin_size)
	bin2 = int(int(e)/bin_size)
	bin_code = str(bin1 * bin_size) + '-' + str((bin1+1) * bin_size)
	if bin1 == bin2: return [bin_code]
	else: return [ bin_code, str(bin2 * bin_size) + '-' + str((bin2+1) * bin_size) ]

def is_contained_in(s, e, x, y):

	s, e = natsorted((s, e))
	x, y = natsorted((x, y))
	return all([s <= x <= e, s <= y <= e])

def extract_ref_features(genome_gff, repr_file, output_dir="./", GFF_feature="CDS"):

	# Reference transcripts
	repr_gene_models = set([ line.rstrip().split('\t')[0] for line in open(repr_file) ])

	# Extract features from reference annotation
	trs, exon2id, bins = set(), {}, { '+': {}, '-': {} }
	with open(output_dir+'repr_CDS.bed', 'w') as fout:
		for line in open(genome_gff):
			if not line.startswith('#'):
				c, source, feature, s, e, score, strand, phase, attr = line.rstrip().split('\t')
				s, e = int(s), int(e)
				if GFF_feature == feature:
					try:

						t_id = re.search('transcript_id "(.*?)";', attr).group(1)
						if t_id in repr_gene_models:
							if t_id in trs: feature_count += 1
							else:
								trs.add(t_id)
								feature_count = 1

							try: exon2id['{}:{}-{}'.format(c, s, e)].append('{}-{}'.format(t_id, feature_count))
							except KeyError: exon2id['{}:{}-{}'.format(c, s, e)] = ['{}-{}'.format(t_id, feature_count)]

							if c not in bins[strand]: bins[strand][c] = {}
							bin_codes = get_bin(s, e, 5000)
							for bin_code in bin_codes:
								try: bins[strand][c][bin_code].add((s, e))
								except KeyError: bins[strand][c][bin_code] = set([(s, e)])

							fout.write( '{}\t{}\t{}\t{}-{}\t1000\t{}\t{}\n'.format(c, s, e, t_id, feature_count, strand, phase) )

					except AttributeError: pass

	return trs, exon2id, bins

def prepare_junction_lib(genome_gff, repr_file, junction_files, sample_names, output_dir="./", exitron_bed=None, min_support=3, GFF_feature="CDS"):

	# Check the number of files
	if len(junction_files) != len(sample_names):
		sys.stderr.write("The number of junction files and sample names should be the same! Aborting...")
		sys.exit(1)

	# Get the shared sample junction
	junction_sets = []
	uniq_names = set(sample_names)
	for n in natsorted(uniq_names):
		indices = [ i for i, x in enumerate(sample_names) if x == n ]
		junction_sets.append( get_shared_junctions([ junction_files[i] for i in indices ]) )
	all_junctions = set.union(*junction_sets)

	# Add existing exitrons when supplied
	if os.path.isfile(exitron_bed):
		existing_EI = set([ ( r['c'], r['s'], r['e']-1, r['strand'] ) for r in yield_bed(exitron_bed) ])
		all_junctions = all_junctions.union(existing_EI)

	# Extract features from reference annotation
	trs, exon2id, bins = extract_ref_features(genome_gff, repr_file, output_dir, GFF_feature)

	# Map junctions to features
	for j in natsorted(all_junctions):
		
		c, s, e, strand = j
		candidates, match = {}, None
		bin_codes = get_bin(s, e, 5000)

		for bin_code in bin_codes:
			try:
				if bin_code in bins[strand][c]:
					for x, y in natsorted(bins[strand][c][bin_code]):

						if is_contained_in(x, y, s, e):
							exon = exon2id['{}:{}-{}'.format(c, x, y)][0]
							candidates[exon] = { 'start': int(x), 'end': int(y) }
			except KeyError: pass

		if len(candidates) == 1: match = ''.join(candidates.keys())
		elif len(candidates) > 1: 
			feature_lengths = { exon: candidates[exon]['end'] - candidates[exon]['start'] for exon in candidates }
			match = max(feature_lengths, key=feature_lengths.get)

		if match != None:
			print '{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(c, candidates[match]['start'], candidates[match]['end'], strand, match, s, e)

def get_EI_from_introns(genome_gff, repr_file, introns_bed, output_dir="./", GFF_feature="CDS"):

	# Extract features from reference annotation
	trs, exon2id, bins = extract_ref_features(genome_gff, repr_file, output_dir, GFF_feature)

	# Map introns to exons
	track = set()
	with open(output_dir+'AtRTD2_intron_map.txt', 'w') as fout:
		for line in open(introns_bed):
			c, s, e, t_exon, score, strand = line.rstrip().split('\t')
			s, e = int(s), int(e)
			candidates, match = {}, None
			bin_codes = get_bin(s, e)

			for bin_code in bin_codes:
				try: 
					if bin_code in bins[strand][c]:
						for x, y in natsorted(bins[strand][c][bin_code]):
							if is_contained_in(x, y, s, e):
								exon = exon2id['{}:{}-{}'.format(c, x, y)][0]
								candidates[exon] = { 'start': int(x), 'end': int(y) }
				except KeyError: pass

			if len(candidates) == 1: match = ''.join(candidates.keys())
			elif len(candidates) > 1: 
				feature_lengths = { exon: candidates[exon]['end'] - candidates[exon]['start'] for exon in candidates }
				match = max(feature_lengths, key=feature_lengths.get)

			if match != None and (s, e) not in track:
				fout.write( '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(c, candidates[match]['start'], candidates[match]['end'], strand, match, s, e) )
				track.add((s, e))

def get_ASta_retained_introns(genome_gff, repr_file, asta_output, output_dir="./", GFF_feature="CDS"):

	# Extract features from reference annotation
	trs, exon2id, bins = extract_ref_features(genome_gff, repr_file, output_dir, GFF_feature)

	# Get retained introns from AStalavista output
	IR_pattern = re.compile('0,\d\^\d-')
	MX_IR_pattern = re.compile('1\^2-,3\^4')
	IR_coor_pattern = re.compile('(\d+)\^(\d+)')
	A3_IR_pattern = re.compile('1-,2-3\^4-')
	A5_IR_pattern = re.compile('1\^2-4\^,3\^$')

	track = set()

	with open(output_dir+'AStalavista_IR_map.txt', 'w') as fout:
		for line in open(asta_output):
			c, source, feature, s, e, score, strand, frame, attributes = line.rstrip().split('\t')
			attr = {}
			for a in attributes.split(';'):
				if len(a):
					attr_name, attr_value = a.split(' "')
					attr[attr_name.strip()] = attr_value.replace('\"', '')

			# Only select IR events
			if any([ IR_pattern.match(attr['structure']),
					MX_IR_pattern.match(attr['structure']),
					A3_IR_pattern.match(attr['structure']),
					A5_IR_pattern.match(attr['structure'])
			]):

				# Iterate the potential multiple IRs
				for (s, e) in re.findall(IR_coor_pattern, attr['splice_chain']):

					s, e = natsorted((int(s), int(e)))
					e = e-1
					candidates, match = {}, None
					bin_codes = get_bin(s, e)

					for bin_code in bin_codes:
						try: 
							if bin_code in bins[strand][c]:
								for x, y in natsorted(bins[strand][c][bin_code]):
									if is_contained_in(x, y, s, e):
										exon = exon2id['{}:{}-{}'.format(c, x, y)][0]
										candidates[exon] = { 'start': int(x), 'end': int(y) }
						except KeyError: pass

				if len(candidates) == 1: match = ''.join(candidates.keys())
				elif len(candidates) > 1:
					feature_lengths = { exon: candidates[exon]['end'] - candidates[exon]['start'] for exon in candidates }
					match = max(feature_lengths, key=feature_lengths.get)

				if match != None and (s, e) not in track:
					fout.write( '{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(c, candidates[match]['start'], candidates[match]['end'], strand, match, s, e) )
					track.add((s, e))

def filter_exitrons(introns_bed, exons_bed, cds_bed, junction_map):
	
	# Get annotated introns
	annotated_introns = set([ '{}:{}-{}'.format(d['c'], d['s'], d['e']) for d in yield_bed(introns_bed) ])

	# Get reference exons
	exons = {}
	for d in yield_bed(exons_bed):
		t_id = d['id'].split('-')[0]
		try: exons[t_id].append([ d['c'], d['s'], d['e'] ])
		except KeyError: exons[t_id] = [[ d['c'], d['s'], d['e'] ]]
	for t_id in exons: exons[t_id].sort(key=lambda x: x[1]) # Sort coordinates

	# Get representative gene models CDS
	repr_CDS = {}
	for d in yield_bed(cds_bed):
		t_id = d['id'].split('-')[0]
		try: repr_CDS[t_id].append([ d['c'], d['s'], d['e'], d['strand'], d['phase'] ])
		except KeyError: repr_CDS[t_id] = [[ d['c'], d['s'], d['e'], d['strand'], d['phase'] ]]
	for t_id in repr_CDS: repr_CDS[t_id].sort(key=lambda x: x[1]) # Sort coordinates

	for line in open(junction_map):
		c, s, e, strand, t_exon, j_start, j_end = line.rstrip().split('\t')
		IEx3 = 'yes' if abs(int(j_end)-int(j_start)) % 3 == 0 else 'no'
		t_id, exon_number = t_exon.split('-')
		g_id = t_id.split('.')[0]

		if '{}:{}-{}'.format(c, j_start, j_end) in annotated_introns:
			
			# Get translation info
			t_id = t_exon.split('-')[0]
			strand = repr_CDS[t_id][0][3]
			if strand == '+':
				phase = repr_CDS[t_id][0][4]
				translation_start = repr_CDS[t_id][0][1]
				first_exon_start = 0
			else:
				phase = repr_CDS[t_id][-1][4]
				translation_start = repr_CDS[t_id][-1][2]

			# Generate exitron spliced in bed file
			with open("tmp_spliced_in.bed", 'w') as fout:
				for c, s, e in exons[t_id]:
					if strand == '+' and not translation_start > e:
						if first_exon_start == 0: first_exon_start = s
						fout.write( '{}\t{}\t{}\n'.format(c, s, e) )
					elif strand == '-' and not translation_start < s:
						first_exon_start = e
						fout.write( '{}\t{}\t{}\n'.format(c, s, e) )

			# Generate exitron spliced out bed file
			with open("tmp_spliced_out.bed", 'w') as fout:
				for c, s, e in exons[t_id]:
					if strand == '+' and translation_start < e or strand == '-' and translation_start > s:
						if s <= int(j_start) <= e and s <= int(j_end) <= e:
							fout.write( '{}\t{}\t{}\n'.format(c, s, j_start) )
							fout.write( '{}\t{}\t{}\n'.format(c, j_end, e) )
						else: fout.write( '{}\t{}\t{}\n'.format(c, s, e) )

			# Get sequences
			subprocess.call("bedtools getfasta -fi $GENOME_FASTA -bed tmp_spliced_in.bed -fo tmp_spliced_in.fa", shell=True)
			subprocess.call("bedtools getfasta -fi $GENOME_FASTA -bed tmp_spliced_out.bed -fo tmp_spliced_out.fa", shell=True)
			i = ''.join([ record.seq for record in yield_fasta('tmp_spliced_in.fa') ])
			o = ''.join([ record.seq for record in yield_fasta('tmp_spliced_out.fa') ])

			# Get translation offset and get revcom if required
			if strand == '+': translation_offset = int(translation_start) - int(first_exon_start) -1
			elif strand == '-':
				i = Sequence(i).get_reverse_complement()
				o = Sequence(o).get_reverse_complement()
				translation_offset = int(first_exon_start) - int(translation_start)
	
			# Get protein sequence
			inProtein = Sequence(i[translation_offset:]).translate()
			try: inProtein = inProtein[:inProtein.index('*')]
			except ValueError: inProtein = inProtein[::]

			outProtein = Sequence(o[translation_offset:]).translate()
			try: outProtein = outProtein[:outProtein.index('*')]
			except ValueError: outProtein = outProtein[::]

			# Check if exitron containing protein (i.e. spliced in variant)
			# is at least longer than the spliced out variant and that the 
			# stop codon is at the same or an later position as with the 
			# spliced out variant.
			keep = False
			if len(inProtein) >= len(outProtein):
				if len(inProtein) != len(outProtein):
					if len(inProtein) >= len(outProtein) + abs(int(j_end)-int(j_start)-1)/3: keep = True
				else: keep = True

			if keep:
				print '{}\tyes\t{}'.format(line.rstrip(), IEx3)

		else: 
			print '{}\tno\t{}'.format(line.rstrip(), IEx3)

	# Clean up
	subprocess.call("rm tmp_spliced_*", shell=True)

def prepare_bam_files(bam_file, output_dir, handle):
	''' Extract the unique reads and unique exonic reads from the supplied bam
	file, index and output to output_dir. '''

	# Create output_dir if not exists
	if not os.path.exists(output_dir): os.makedirs(output_dir)

	# Check if files already exists
	a1 = True
	uniq_reads_bam = '{}accepted_hits_uniq_reads.{}.bam'.format(output_dir, handle)
	if os.path.exists(uniq_reads_bam):
		a1 = query_yes_no('{} already exists. Do you want to overwrite it?'.format(uniq_reads_bam))

	if a1: 
		cmd = "samtools view %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]=~/N/){ print \"$_\"; }' > %s" % (bam_file, uniq_reads_bam.replace('.bam', '.sam'))
		subprocess.call(cmd, shell=True)

		# Convert sam to bam
		subprocess.call("samtools view -bT $GENOME_FASTA {} > {}".format(uniq_reads_bam.replace('.bam', '.sam'), uniq_reads_bam), shell=True)

		# Index bam
		subprocess.call("samtools index {}".format(uniq_reads_bam), shell=True)

		# Remove sam
		subprocess.call("rm %s" % (uniq_reads_bam.replace('.bam', '.sam')), shell=True)

	a2 = True
	uniq_exon_reads_bam = '{}accepted_hits_uniq_exonic_reads.{}.bam'.format(output_dir, handle)
	if os.path.exists(uniq_exon_reads_bam):
		a2 = query_yes_no('{} already exists. Do you want to overwrite it?'.format(uniq_exon_reads_bam))

	if a2:
		cmd = "samtools view %s | grep -w \"NH:i:1\" | perl -n -e '@line=split(/\\t/,$_); if ($line[5]!~/N/){ print \"$_\"; }' > %s" % (bam_file, uniq_exon_reads_bam.replace('.bam', '.sam'))
		subprocess.call(cmd, shell=True)

		# Convert sam to bam
		subprocess.call("samtools view -bT $GENOME_FASTA {} > {}".format(uniq_exon_reads_bam.replace('.bam', '.sam'), uniq_exon_reads_bam), shell=True)

		# Index bam
		subprocess.call("samtools index {}".format(uniq_exon_reads_bam), shell=True)

		# Remove sam
		subprocess.call("rm %s" % (uniq_exon_reads_bam.replace('.bam', '.sam')), shell=True)

def calculate_PSI(exitron_map, uniq_reads, uniq_exonic):
	''' Calculate exitron PSI values, based on the coverage of the 
	unique exonic reads. 
						
						A		   B		 C
					   ----		  ----		----
	EEEEEEEEEEEEEEEEEEEEEXXXXXXXXXXXXXXXXXXXXXEEEEEEEEEEEEEEEEEEEEE
					   --	                  --
						 \                   / 
						  \        D        / 							
                           -----------------
                           	
    E = exon, X = exitron
    A = Reads aligning from -10 to +10 around exitron 5'SS
    B = Reads aligning from -10 to +10 around exitron middle point
    C = Reads aligning from -10 to +10 around exitron 3'SS
    D = Reads supporting exitron splicing, output from get_junction_uniq_read_support()
    
    PSI = ( ((A+B+C)/3) / ((A+B+C)/3) + D ) * 100 '''

	# Check if all files are there
	if all([ os.path.exists(exitron_map), os.path.exists(uniq_reads), os.path.exists(uniq_exonic) ]):
		
		# Get junction support
		rc, info = {}, {}
		for line in open(exitron_map):

			c, s, e, strand, ID, j_start, j_end, is_annotated, IEx3 = line.rstrip().split('\t')
			info['{}:{}-{}'.format(c, j_start, int(j_end)+1)] = { 'strand': strand, 'ID': ID, 'is_annotated': is_annotated, 'IEx3': IEx3 }
			subprocess.call("samtools view {} {}:{}-{} > tmp.sam".format(uniq_reads, c, s, e), shell=True)

			# Get the required read/junction gap signature
			N = "%dN" % (int(j_end)-int(j_start))

			uniq_count = 0
			for aln in open("tmp.sam"):
				qname = aln.split('\t')[0]
				pos = int(aln.split('\t')[3])
				cigar = aln.split('\t')[5]
				start = (pos + int(re.search('^([0-9]+)M', cigar).group(1))) - 1
				
				# Check if the junctions is at the correct position 
				# and if the junction size is correct
				if N in cigar and start == int(j_start): uniq_count += 1

			rc['{}:{}-{}'.format(c, j_start, int(j_end)+1)] = { 'A': 0, 'B': 0, 'C': 0, 'D': uniq_count }

		# Get A, B, C read support
		with open("tmp_coverageBed_input.bed", 'w') as fout: 
			for line in open(exitron_map):

				c, s, e, strand, ID, j_start, j_end, is_annotated, IEx3 = line.rstrip().split('\t')
				j_start, j_end = int(j_start), int(j_end)
				middle_point = int(j_start + ((j_end-j_start)/2))
				locus = '{}:{}-{}'.format(c, j_start, int(j_end)+1)

				fout.write( '{}\t{}\t{}\t{}_A\n'.format(c, j_start-10, j_start+10, locus) )
				fout.write( '{}\t{}\t{}\t{}_B\n'.format(c, middle_point-10, middle_point+10, locus) )
				fout.write( '{}\t{}\t{}\t{}_C\n'.format(c, j_end-10, j_end+10, locus) )

		subprocess.call("sort -k1,1 -k2,2n tmp_coverageBed_input.bed > tmp_coverageBed_input.sorted.bed", shell=True)
		subprocess.call("coverageBed -sorted -counts -a tmp_coverageBed_input.sorted.bed -b {} > tmp_coverageBed_output.bed".format(uniq_exonic), shell=True)
		for line in open("tmp_coverageBed_output.bed"):
			c, start, end, locus, coverage = line.rstrip().split('\t')
			locus, letter = locus.split('_')
			rc[locus][letter] = int(coverage)

		# Calculate PSI
		for x in natsorted(rc):
			try: 
				PSI = ( ( ( rc[x]['A'] + rc[x]['B'] + rc[x]['C'] ) / 3 ) / ( ( ( rc[x]['A'] + rc[x]['B'] + rc[x]['C'] ) / 3 ) + rc[x]['D'] ) ) * 100
				print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}'.format(x, info[x]['ID'], info[x]['strand'], info[x]['is_annotated'], info[x]['IEx3'], rc[x]['A'], rc[x]['B'], rc[x]['C'], rc[x]['D'], PSI)
			except ZeroDivisionError: 
				print '{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\tNA'.format(x, info[x]['ID'], info[x]['strand'], info[x]['is_annotated'], info[x]['IEx3'], rc[x]['A'], rc[x]['B'], rc[x]['C'], rc[x]['D'])

		# Clean up
		subprocess.call("rm tmp.sam tmp_coverageBed_*.bed", shell=True)

	else:
		sys.stderr.write("One (or more) input files cannot be accessed.\n")

def parse_exitron_file(fname):

	e_dict = {}
	for line in open(fname):
		locus, t_exon, strand, is_annotated, IEx3, A, B, C, D, PSI = line.rstrip().split('\t')
		t_id, exon_number = t_exon.split('-')
		g_id = t_id.split('.')[0]
		e_dict[locus] = { 
			'g_id': g_id,
			't_id': t_id,
			'exon_number': exon_number,
			'strand': strand,
			'is_annotated': is_annotated,
			'IEx3': IEx3,
			'A': int(A),
			'B': int(B),
			'C': int(C),
			'D': int(D),
			'PSI': float(PSI) if PSI != "NA" else PSI,
			'line': line.rstrip()
		}

	return e_dict

def compare_exitrons(s1_files, s2_files, mode="delta", min_delta=15, max_std=5):

	import numpy as np
	from scipy.stats import ttest_ind

	print 'GENE\tCOORD\tEVENT_SIZE\tEIx3\tAS_TYPE\t{}\t{}\tDIRECTION\tDELTA'.format( '\t'.join([ f.split('/')[-1].split('.')[0] for f in s1_files ]), '\t'.join([ f.split('/')[-1].split('.')[0] for f in s2_files ]) )

	s1, s2 = [ [ parse_exitron_file(f) for f in s ] for s in [s1_files, s2_files] ]
	for locus in natsorted(s1[0]):

		s1_val, s2_val = [ [s[i][locus]['PSI'] for i in xrange(len(s))] for s in [s1, s2] ]
		if not 'NA' in s1_val + s2_val:
			s1_std, s2_std = np.std(s1_val), np.std(s2_val)
			s1_mean, s2_mean = np.mean(s1_val), np.mean(s2_val)
			delta = abs(s1_mean - s2_mean)

			exceeds_expectations = False
			if mode == "delta":
				if all([ delta >= min_delta, s1_std < max_std, s2_std < max_std ]): exceeds_expectations = True

			elif mode == "ttest":
				pval = ttest_ind(s1_val, s2_val, equal_var=False)[1]
				if all([ pval < 0.05, s1_std < max_std, s2_std < max_std, delta >= min_delta ]): exceeds_expectations = True

			if exceeds_expectations:
				c, s, e = re.search('(\d):(\d+)-(\d+)', locus).groups()
				info = s1[0][locus]

				event_size = abs(int(e)-int(s))-1
				print '{}\t{}:{}-{}\t{}\t{}\tEI\t{}\t{}\t{}\t{:.2f}'.format(info['g_id'], c, s, e, event_size, "yes" if event_size % 3 == 0 else "no", '\t'.join(map(str, s1_val)), '\t'.join(map(str, s2_val)), 'CTRL<TEST' if s1_mean < s2_mean else 'CTRL>TEST', delta)

				'''
				print '{}\tlaBarta\texitron\t{}\t{}\t.\t{}\t.\tgene_id "{}"; transcript_id "{}"; exon_number "{}"; is_annotated "{}"; IEx3 "{}"; {} "{}"; PSI_1 "{:.3f}"; PSI_2 "{:.3f}";'.format(
					c, s, e,
					info['strand'],
					info['g_id'],
					info['t_id'],
					info['exon_number'],
					info['is_annotated'],
					info['IEx3'],
					'delta' if mode == "delta" else "pvalue",
					delta if mode == "delta" else pval,
					s1_mean,
					s2_mean
				)
				'''

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	subparsers = parser.add_subparsers(dest='command', help='sub-command help')

	parser_a = subparsers.add_parser('prepare-junctions', help="Select all unique junctions with a minimum read support (default = 3) and map to the CDS (default) of the representative gene models.")
	parser_a.add_argument('-g', '--genome-gff', required=True, help="Genome GTF/GFF annotation file.")
	parser_a.add_argument('-r', '--ref', required=True, help="List of representative gene models (transcripts).")
	parser_a.add_argument('-j', '--junctions', required=True, nargs='+', help="Junction bed files to use.")
	parser_a.add_argument('-n', '--names', required=True, nargs='+', help="List of condition names, e.g. WT WT WT OX OX OX")
	parser_a.add_argument('-o', '--output-dir', default=".", help="Output directory.")
	parser_a.add_argument('-e', '--exitrons', default=None, help="Existing exitron bed file.")
	parser_a.add_argument('-m', '--min-support', type=int, default=3, help="Minimum required junction support.")
	parser_a.add_argument('-f', '--feature', default="CDS", help="GTF/GFF feature to select.")

	parser_b = subparsers.add_parser('AStalavista-EI', help="Select all intron retention events from the AStalavista output and map to the CDS of the representative gene models.")
	parser_b.add_argument('-g', '--genome-gff', required=True, help="Genome GTF/GFF annotation file.")
	parser_b.add_argument('-r', '--ref', required=True, help="List of representative gene models (transcripts).")
	parser_b.add_argument('-a', '--asta', required=True, help="AStalavista stand-alone version asta GTF output file")
	parser_b.add_argument('-o', '--output-dir', default=".", help="Output directory.")
	parser_b.add_argument('-f', '--feature', default="CDS", help="GTF/GFF feature to select.")

	parser_i = subparsers.add_parser('intron-EI', help="Map all introns to the CDS of the representative gene models.")
	parser_i.add_argument('-g', '--genome-gff', required=True, help="Genome GTF/GFF annotation file.")
	parser_i.add_argument('-r', '--ref', required=True, help="List of representative gene models (transcripts).")
	parser_i.add_argument('-i', '--introns', required=True, help="Introns bed file.")
	parser_i.add_argument('-o', '--output-dir', default=".", help="Output directory.")
	parser_i.add_argument('-f', '--feature', default="CDS", help="GTF/GFF feature to select.")

	# Filter exitrons subparser
	parser_c = subparsers.add_parser('filter-exitrons', help="Run after prepare-junctions. Filter out intron retention events from the potential exitrons. Output is redirected to STDOUT.")
	parser_c.add_argument('-i', '--introns', required=True, help="Genome intron bed file.")
	parser_c.add_argument('-e', '--exons', required=True, help="Genome exon bed file.")
	parser_c.add_argument('-c', '--CDS', required=True, help="CDS bed of representative gene models.")
	parser_c.add_argument('-j', '--junction-map', required=True, help="Junction map file (from prepare-junctions).")

	# Prepare bam file subparser
	parser_d = subparsers.add_parser('prepare-bam', help='Extract the unique (exonic) reads from the TopHat2 mapping bam file.')
	parser_d.add_argument('-b', '--bam', required=True, help='TopHat2 bam mapping file.')
	parser_d.add_argument('-o', '--output-dir', required=True, help="Output directory.")
	parser_d.add_argument('-f', '--file-handle', required=True, help="Unique file handle. The files will be accepted_hits_uniq_reads.[file-handle].bam and accepted_hits_uniq_exonic_reads.[file-handle].bam.")

	# Calculate PSI subparser
	parser_e = subparsers.add_parser('calculate-PSI', help="Calculate the PSI. Output is directed to the STDOUT.")
	parser_e.add_argument('--exitron-map', required=True, help='Junctions to coding exon mapping file.')
	parser_e.add_argument('--uniq-reads', required=True, help='Unique reads bam file.')
	parser_e.add_argument('--uniq-exonic', required=True, help='Unique exonic reads bam file.')

	# Compare exitrons subparser
	parser_f = subparsers.add_parser('compare', help="Compare exitron PSIs between conditions. Output is directed to the STDOUT.")
	parser_f.add_argument('-s1', required=True, nargs='+', help="Condition one exitron file regex.")
	parser_f.add_argument('-s2', required=True, nargs='+', help="Condition two exitron file regex.")
	parser_f.add_argument('--mode', choices=["delta", "ttest"], default="delta", help="Mode to compare exitrons.")
	parser_f.add_argument('--min-delta', default=15, type=int, help="Minimal mean difference between samples.")
	parser_f.add_argument('--max-std', default=5, type=int, help="Maximal internal sample standard deviation")

	args = parser.parse_args()
	if args.command == "prepare-junctions":
		prepare_junction_lib(args.genome_gff, args.ref, args.junctions, args.names, add_slash(args.output_dir), args.exitrons, args.min_support, args.feature)
	elif args.command == "AStalavista-EI":
		get_ASta_retained_introns(args.genome_gff, args.ref, args.asta, add_slash(args.output_dir), args.feature)
	elif args.command == "intron-EI":
		get_EI_from_introns(args.genome_gff, args.ref, args.introns, add_slash(args.output_dir), args.feature)
	elif args.command == "filter-exitrons":
		filter_exitrons(args.introns, args.exons, args.CDS, args.junction_map)
	elif args.command == "prepare-bam":
		prepare_bam_files(args.bam, args.output_dir, args.file_handle)
	elif args.command == "calculate-PSI":
		calculate_PSI(args.exitron_map, args.uniq_reads, args.uniq_exonic)
	elif args.command == "compare":
		compare_exitrons(args.s1, args.s2, args.mode, args.min_delta, args.max_std)