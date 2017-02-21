#!/usr/bin/env python

"""
Collection of often used fileparser
"""

def yield_salmon(f):
	""" Yield salmon output per line """

	for line in open(f):
		try: 
			t_id, length, eff_length, tpm, num_reads = line.rstrip().split('\t')
			yield( { 'transcript_id': t_id,  \
				'length': int(length), \
				'eff_length': float(eff_length), \
				'tpm': float(tpm), \
				'num_reads': float(num_reads) } )
		except ValueError: pass

def parse_salmon(f, group="gene"):
	""" Return the TPM values and group per gene or transcript """

	if group not in ["gene", "transcript"]: group = "gene"
	
	salmon = {}
	for record in yield_salmon(f):
		t_id = record["transcript_id"]
		if group == "gene":
			g_id = t_id.split('.')[0]
			try: salmon[g_id][t_id] = float(record["tpm"])
			except KeyError: salmon[g_id] = { t_id: float(record["tpm"]) }
		else: 
			salmon[t_id] = float(record["tpm"])

	return salmon

def yield_kallisto(f):
	""" Yield kallisto output per line """

	for line in open(f):
		try: 
			t_id, length, eff_length, est_counts, tpm = line.rstrip().split('\t')
			yield( { 'transcript_id': t_id,  \
				'length': int(length), \
				'eff_length': float(eff_length), \
				'tpm': float(tpm), \
				'est_counts': float(num_reads) } )
		except ValueError: pass

	pass

def parse_kallisto(f, group="gene"):
	""" Return the TPM values and group per gene or transcript """

	if group not in ["gene", "transcript"]: group = "gene"
	
	kallisto = {}
	for record in yield_kallisto(f):
		t_id = record["transcript_id"]
		if group == "gene":
			g_id = t_id.split('.')[0]
			try: kallisto[g_id][t_id] = float(record["tpm"])
			except KeyError: kallisto[g_id] = { t_id: float(record["tpm"]) }
		else: 
			kallisto[t_id] = float(record["tpm"])

	return kallisto

def yield_fasta(f):
	''' Simple fasta parser that yield's the identifier and sequence of each record '''
	
	class SeqRecord:
		def __init__(self, seq_id, seq):
			self.id = seq_id
			self.seq = seq

	seq = ''
	for line in open(fasta_file):
		if line.startswith('>'):
			if len(seq):
				yield(SeqRecord(identifier, seq))
				seq = ''
			identifier = line[1:].rstrip()
		else: seq += line.rstrip()
	yield(SeqRecord(identifier, seq))