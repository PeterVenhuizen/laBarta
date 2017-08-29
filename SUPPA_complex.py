#!/usr/env/python

"""
SUPPA_complex.py

Goal 1: Identify all the simple and complex events annotated by SUPPA. 
Goal 2: Filter on expressed transcripts
"""

from natsort import natsorted
import argparse

def parse_salmon():
	pass

class Event(object):

	def __init__(self, event_id, impact_region, inc, exc):
		self.event_id = event_id
		self.gene_id = event_id.split(';')[0]
		self.event_type = event_id.split(';')[1].split(':')[0]
		self.impact_region = impact_region
		self.inc = inc
		self.exc = exc

	def is_expressed(self, ctrl_TPM, smpl_TPM, N):
		
		all_true = []
		for i in xrange(0, N):
			if sum([ ctrl_TPM[x][i] for x in self.inc ]) > 1 or sum([ ctrl_TPM[x][i] for x in self.exc ]) > 1 and sum([ smpl_TPM[x][i] for x in self.inc ]) > 1 or sum([ smpl_TPM[x][i] for x in self.exc ]) > 1:
				all_true.append(True)
			else: all_true.append(False)

		return all(all_true)

def yield_SUPPA(ioe_file, boundary_type="S"):
	
	for line in open(ioe_file):
		seqname, gene_id, event_id, alt_transcripts, total_transcripts = line.rstrip().split('\t')

		try: 
			event_type = event_id.split(';')[1].split(':')[0]

			if event_type == "A3":
				event_type, c, short_exon, long_exon, strand = event_id.split(';')[1].split(':')
				if strand == '+':
					e1, s2 = map(int, short_exon.split('-'))
					e1, s3 = map(int, long_exon.split('-'))

				elif strand == '-':
					e2, s3 = map(int, short_exon.split('-'))
					e1, s3 = map(int, long_exon.split('-'))
				impact_region = [e1, s3]

			elif event_type == "A5":
				event_type, c, short_exon, long_exon, strand = event_id.split(';')[1].split(':')
				if strand == '+':
					e2, s3 = map(int, short_exon.split('-'))
					e1, s3 = map(int, long_exon.split('-'))

				elif strand == '-':
					e1, s2 = map(int, short_exon.split('-'))
					e1, s3 = map(int, long_exon.split('-'))
				impact_region = [e1, s3]

			elif event_type == "MX":
				pass

			elif event_type == "SE":
				event_type, c, intron_one, intron_two, strand = event_id.split(';')[1].split(':')
				e1, s2 = map(int, intron_one.split('-'))
				e2, s3 = map(int, intron_two.split('-'))
				impact_region = [e1, s3]

			elif event_type in ["RI", "EI"]:
				event_type, c, exon_left, intron, exon_right, strand = event_id.split(';')[1].split(':')
				impact_region = map(int, intron.split('-'))

			# Yield event
			try: 
				yield(Event(event_id, impact_region, alt_transcripts.split(','), total_transcripts.split(',')))
			except UnboundLocalError: pass

		except IndexError: pass # Skip header

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

def SUPPA_complex(ioe_files, ctrl_TPM_file, smpl_TPM_file, boundary_type="S"):

	# Parse TPMs
	ctrl_TPM = {}
	for line in open(ctrl_TPM_file):
		try: ctrl_TPM[line.split('\t')[0]] = map(float, line.rstrip().split('\t')[1:])
		except ValueError: N = len(line.rstrip().split('\t'))

	smpl_TPM = {}
	for line in open(smpl_TPM_file):
		try: smpl_TPM[line.split('\t')[0]] = map(float, line.rstrip().split('\t')[1:])
		except ValueError: pass

	# Group SUPPA events per gene
	events = {}
	for f in ioe_files:
		for e in yield_SUPPA(f):

			#if e.is_expressed(ctrl_TPM, smpl_TPM, N):

			try: events[e.gene_id].append(e)
			except KeyError: events[e.gene_id] = [e]

	# Loop through the genes
	for g_id in natsorted(events):

		# Is there more than 1 event?
		if len(events[g_id]) > 1:
			
			# Yes -> Look for overlapping events

			# Compare the first event against the rest
			gene_events = set([ e for e in events[g_id] ])
			
			while len(gene_events) > 0:
				magnitude = gene_events.pop() # Pop! Pop!				
				pop_pop = set([ e for e in gene_events if is_overlapping(magnitude.impact_region[0], magnitude.impact_region[1], e.impact_region[0], e.impact_region[1]) ])

				# Are two or more events overlapping?
				if len(pop_pop) == 0: # No -> Simple
					print '{}\tSIMPLE\t{}'.format(g_id, magnitude.event_id)

				elif len(pop_pop) == 1: # Might be strict IR event with the same intron

					pop = next(iter(pop_pop))
					if pop.impact_region == magnitude.impact_region:
						print '{}\tSIMPLE\t{}'.format(g_id, magnitude.event_id)
						print '{}\tSIMPLE\t{}'.format(g_id, pop.event_id)
					else:
						print '{}\tCOMPLEX\t{},{}'.format(g_id, magnitude.event_id, pop.event_id)

				else:
					print '{}\tCOMPLEX\t{},{}'.format(g_id, magnitude.event_id, ','.join([ e.event_id for e in pop_pop ]))

				for e in pop_pop: gene_events.remove(e)

				# Yes -> Are the transcripts involved actually expressed?

					# Filter on expressed transcripts and filter the 
					# SUPPA events accordingly 

					# Are two or more events overlapping?

						# Yes -> Complex

						# No -> Simple

		# No -> It must be a simple event
		else:
			print '{}\tSIMPLE\t{}'.format(g_id, events[g_id][0].event_id)

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('--ioe', nargs='+', required=True, help="SUPPA ioe files.")
	parser.add_argument('--control-TPM', required=True, help="Control TPM values.")
	parser.add_argument('--sample-TPM', required=True, help="Sample TPM values.")
	parser.add_argument('--boundary-type', choices=["S", "V"], default="S", help="SUPPA generateEvents boundary type.")
	args = parser.parse_args()

	SUPPA_complex(args.ioe, args.control_TPM, args.sample_TPM, args.boundary_type)