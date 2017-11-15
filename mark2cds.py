#!/usr/bin/env python

"""
mark2cds.py -> Translate Mark's protein translation to CDS annotation
"""

import re
import argparse
from natsort import natsorted
from fileParser import yield_fasta
from sequence import Sequence

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

def mark2cds(fasta_file, exons_fasta):
	
	marks = { record.id: record.seq for record in yield_fasta(fasta_file) }
	exons = parse_fasta_file({}, exons_fasta)

	ORF = re.compile('(M.*?\*)')
	for g_id in natsorted(exons):

		isFwd = exons[g_id]['isFwd']
		del exons[g_id]['isFwd']
		for t_id in natsorted(exons[g_id]):

			if t_id in marks:

				#print t_id
				track_len = 0

				exs = [ e for e in natsorted(exons[g_id][t_id]) ] if isFwd else [ e for e in natsorted(exons[g_id][t_id], reverse=True) ]
				mark = marks[t_id]
				#print mark

				if isFwd: rna = ''.join([ exons[g_id][t_id][e] for e in natsorted(exons[g_id][t_id]) ])
				else: rna = ''.join([ exons[g_id][t_id][e] for e in natsorted(exons[g_id][t_id], reverse=True) ])

				for i in xrange(0, 3):
					trans, leftover = Sequence(rna).translate(i)
					if mark in trans:
						m2star = trans[trans.index(mark):]
						try: m2star = ORF.search(m2star).group(0)
						except AttributeError: pass
						#print m2star, leftover

						start_found, end_found = False, False

						for i, ex in enumerate(exs):
							s, e = map(int, ex.split('-'))
							for j in xrange(0, 3):

								ex_trs, leftover = Sequence(exons[g_id][t_id][ex]).translate(j)
								#print ex, j, ex_trs, leftover

								if ex_trs in m2star:

									# Check if maybe the translation maybe exactly matches the protein start
									if not start_found and m2star.startswith(ex_trs): 
										start_found = True
										s += j

									print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], s, e, '+' if isFwd else '-', j, t_id, g_id)
									track_len += len(ex_trs) + bool(leftover)
									#print ex, j, s, e, ex_trs, leftover
									#print track_len, len(ex_trs)+bool(leftover), ex_trs, leftover
									break

								else:

									if not start_found:
										
										for n in xrange(0, len(ex_trs)):

											# Check whether the start and end of the CDS are in the same exon
											try: m = n + ex_trs[n:].index('*') + 1
											#except ValueError: m = -1
											except ValueError: m = len(ex_trs)
											#print m

											if m2star.startswith(ex_trs[n:m]) and len(ex_trs[n:m]) > 1:

												#print 'START', ex_trs, ex_trs[n:m], leftover

												start_found = True
												#if m != -1: end_found = True
												if m != len(ex_trs): end_found = True

												if start_found and end_found:
													if isFwd: 
														start = s+(n*3)+j
														end = s + (m*3)-1+j
													else:
														end = s+(len(ex_trs[n:])*3)-1+len(leftover)
														#start = end - (len(ex_trs[n:m])*3)
														start = s+(len(ex_trs[m:]*3)-1)+len(leftover)+1
													print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], start, end, '+' if isFwd else '-', j, t_id, g_id)
												elif isFwd: 
													start = s+(n*3)+j
													print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], start, e, '+' if isFwd else '-', j, t_id, g_id)
												else: 
													print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], s, s+(len(ex_trs[n:])*3)-1+len(leftover), '+' if isFwd else '-', j, t_id, g_id)
												track_len += len(ex_trs[n:m]) + bool(leftover)
												#print track_len, len(ex_trs)+bool(leftover), ex_trs, leftover
												break

										# Special case where the start codon 
										# is split over two exons
										if any([not start_found and leftover == 'A', not start_found and leftover == 'AT']):

											#print 'YEAH'
											try: 
												next_ex = exs[i+1]
												for k in xrange(0, 3):
													next_trs, left = Sequence(exons[g_id][t_id][next_ex]).translate(k)
													if next_trs in m2star:
														if not isFwd: s, e = e, s
														if leftover == 'A' and exons[g_id][t_id][next_ex][:k] == "TG":
															print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], e-1, e+1, '+' if isFwd else '-', k, t_id, g_id)
															start_found = True
															track_len += 1
														elif leftover == "AT" and exons[g_id][t_id][next_ex][:k] == "G":
															print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], e, e+1, '+' if isFwd else '-', k, t_id, g_id)
															start_found = True
															track_len += 1
											except IndexError: pass

									elif not end_found:
										
										for n in xrange(len(ex_trs), 0, -1):
											#if m2star.endswith(ex_trs[:n]): print track_seq+ex_trs[:n], leftover
											#if m2star.endswith(ex_trs[:n]) and track_seq+ex_trs[:n] == mark:
											#if m2star.endswith(ex_trs[:n]): print len(mark), track_len+len(ex_trs[:n])
											if m2star.endswith(ex_trs[:n]) and (track_len+len(ex_trs[:n])) == len(mark):

												#print 'END', ex_trs[:n]

												if isFwd: print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], s, s+(n*3)-1+j, '+' if isFwd else '-', j, t_id, g_id)
												else: 

													#start = s+(n*3)+j
													#start = e-(len(ex_trs[:n])*3)
													start = s+(len(ex_trs[n:])*3)+len(leftover)
													print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t{}\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], start, e, '+' if isFwd else '-', j, t_id, g_id)

												#print ex, j, s, s+(n*3)-1+j, ex_trs[:n]
												#print track_len, len(ex_trs)+bool(leftover), ex_trs, leftover
												end_found = True
												break

									#print start_found, end_found

							# Second round for special case start detection
							# NEED TO FIX AT1G03860.P3!!!
							if not start_found:
								for j in xrange(0, 3):
									ex_trs, leftover = Sequence(exons[g_id][t_id][ex]).translate(j)

									for n in xrange(0, len(ex_trs)):

										# Check for the case in which the CDS start
										# is the last AA of an exon
										if m2star.startswith(ex_trs[n:]):
											start_found = True
											if isFwd:
												start = s+(n*3)+j
												print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t+\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], start, end, j, t_id, g_id)
											else:
												print '{}\tMark2CDS\tCDS\t{}\t{}\t.\t-\t{}\ttranscript_id "{}"; gene_id "{}";'.format(t_id[2:3], s, s+(len(ex_trs[n:])*3)-1+len(leftover), j, t_id, g_id)
											break
						break

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--fasta', required=True, help="Mark's translations fasta.")
	parser.add_argument('-e', '--exons', required=True, help="Exons fasta file.")
	args = parser.parse_args()

	mark2cds(args.fasta, args.exons)