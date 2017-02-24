#!/usr/bin/env python

from natsort import natsorted
import sys
sys.path.append('/home/venhuip8/scripts')

def add_slash(path):
	return path if path.endswith('/') else path + '/'

def parse_fasta_file(d, fasta_file):
	
	from fileParser import yield_fasta

	for seq_record in yield_fasta(fasta_file):
		
		t_id = seq_record.id.split('|')[0]
		g_id = t_id.split('.')[0]
		locus = seq_record.id.split(' ')[4].split(':')[1]
		isFwd = True if "FORWARD" in seq_record.id.upper() else False
		
		if g_id not in d: d[g_id] = { 'isFwd': isFwd }
		try: d[g_id][t_id][locus] = seq_record.seq
		except KeyError: d[g_id][t_id] = { locus: seq_record.seq }
	
	return d

def is_overlapping(s, e, x, y):
	''' Check if s-e overlaps with x-y '''

	s, e = natsorted((s, e))
	x, y = natsorted((x, y))

	return any([
		s <= x <= e and s <= y <= e, # Ref exon retained in trs exon
		x <= s <= y and x <= e <= y, # Trs exon retain in ref exon
		x <= s <= y and x < e, # Right overlap with ref exon
		s < x and x <= e <= y # Left overlap with ref exon
	])