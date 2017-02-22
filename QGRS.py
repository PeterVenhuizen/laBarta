#!/usr/bin/env python

"""
Predict the presence of Quadruplex forming G-rich sequences (QGRS) 
in nucleotides sequences (http://bioinformatics.ramapo.edu/QGRS/).
"""

import argparse
from fileParser import yield_fasta

def QGRS(seq, MAX_GRS_LENGTH=45, GGROUPSIZE=2, LOOP_MIN=0, LOOP_MAX=36):
	MIN_GRS_LENGTH = 10
	seqlen = len(seq)
	count = 0
	myset = []
	chosen = []

	# goes through all possible G tract
	for trav in xrange(seqlen-MIN_GRS_LENGTH+1):
			
		# keep going until G found
		if seq[trav] != 'G': continue

		# goes through all possible Gquad
		for findlen in xrange(MIN_GRS_LENGTH, MAX_GRS_LENGTH+1):
			
			# goes through all possible G tract
			for findGlen in xrange(GGROUPSIZE, 7):
				
				if findGlen*4+2 > findlen: break
				
				# check if there is GG, GGG, GGGG in the beginning, first g-tract
				flag = 0
				for checkG in xrange(trav, trav+findGlen):
					if seq[checkG] != 'G':
						flag = 1
						break
				if flag == 1: break
				
				# check if there is GG, GGG, GGGG in the end
				flag = 0
				checkG = trav+findlen-1
				while checkG >= trav+findlen-findGlen: # checks for last g-tract
					if checkG >= seqlen or seq[checkG] != 'G':
						flag = 1
						break
					checkG = checkG-1
				if flag == 1: break
				
				for checkPos2 in xrange(findGlen, findlen-findGlen*3): # check for 2nd g-tract
					
					# check if there is GG, GGG, GGGG at pos2
					flag = 0
					for checkG in xrange(checkPos2+trav, trav+findGlen+checkPos2):
						if seq[checkG] != 'G':
							flag = 1
							break
					if flag == 1: continue
					
					# check if loop size 1 is within limits
					loop_size = checkPos2 - findGlen
					loopSize1 = loop_size
					if loop_size < LOOP_MIN or loop_size > LOOP_MAX: continue
					if checkPos2 == findGlen: zero_loop = 1
					else: zero_loop = 0
					
					# if we had a zero loop, make next one non-zero
					start3 = findGlen+checkPos2 if not zero_loop else findGlen+checkPos2+1
					
					# check for GG, GGG, GGGG at pos3		
					for checkPos3 in xrange(start3, (findlen-findGlen*2)+1):
						flag = 0
						for checkG in xrange(trav+checkPos3, trav+findGlen+checkPos3):
							if seq[checkG] != 'G':
								flag = 1
								break
						if flag == 1: continue
						
						# check if loop#2 is within limits
						loop_size = checkPos3-checkPos2-findGlen
						loopSize2 = loop_size
						if loop_size < LOOP_MIN or loop_size > LOOP_MAX: continue
						
						# check if loop#3 is within limits
						loop_size = findlen-checkPos3-findGlen*2
						loopSize3 = loop_size
						if loop_size < LOOP_MIN or loop_size > LOOP_MAX: continue
						
						# zero loop counter update
						if checkPos3 == checkPos2+findGlen: zero_loop = zero_loop + 1
						if checkPos2 == findGlen: zero_loop = zero_loop + 1
						
						# ensure only one zero loop max!
						if zero_loop == 2: break
							
						Gtract1Start = trav+1
						loopStart1 = Gtract1Start+findGlen

						Gtract2Start = loopStart1+loopSize1
						loopStart2 = Gtract2Start+findGlen
						
						Gtract3Start = loopStart2+loopSize2
						loopStart3 = Gtract3Start+findGlen

						Gtract4Start = loopStart3+loopSize3
						
						g1 = checkPos2-findGlen
						g2 = checkPos3-checkPos2-findGlen
						g3 = findlen-findGlen-checkPos3-findGlen
						
						gap1 = abs(g1-g2)
						gap2 = abs(g1-g3)
						gap3 = abs(g3-g2)
						
						# Calculate score
						FORMULA_LENGTH = MAX_GRS_LENGTH-9
						nscore = FORMULA_LENGTH-(gap1+gap2+gap3)/2+(findGlen-2)*FORMULA_LENGTH
						
						if nscore < 0: continue
						elif loopSize1 == 0 or loopSize2 == 0 or loopSize3 == 0: continue
						else:
							
							lineG = []
							for goThrough in xrange(Gtract1Start-1, Gtract4Start+findGlen-1):
								lineG += seq[goThrough]
							
							aa = { 
								'id': count,
								'glen': findGlen,
								'len': findlen,
								'chosen': 0,
								'totalCount': 0,
								'nscore': nscore,
								'pos': Gtract1Start,
								'checkPos2': checkPos2,
								'checkPos3': checkPos3,
								'g1Start': Gtract1Start,
								'g2Start': Gtract2Start,
								'g3Start': Gtract3Start,
								'g4Start': Gtract4Start,
								'l1size': loopSize1,
								'l1start': loopStart1,
								'l2size': loopSize2,
								'l2start': loopStart2,
								'l3size': loopSize3,
								'l3start': loopStart3,
								'gEnd': Gtract4Start+findGlen-1,
								'GQuadSeq': ''.join(lineG)
							}
							
							count = count+1
							myset.append(aa)
					
					if checkPos3 == checkPos3+findlen:
						zero_loop = zero_loop-1

	# Sort on highest nscore in descending order
	myset = sorted(myset, key=lambda k: k['nscore'], reverse=True)

	# Create masking list
	mask = [0]*seqlen

	# Find all non-overlapping QGRS
	chosencnt = 0
	for m in xrange(len(myset)):
		
		# Check if QGRS is unique
		flag = 0
		for n in xrange(myset[m]['pos'], myset[m]['pos']+myset[m]['len']):
			try: 
				if mask[n] == 1:
					flag = 1
					break
			except IndexError: pass

		# If unique, mask sequence
		if flag == 0:
			chosencnt = chosencnt+1
			myset[m]['chosen'] = 1
			for n in xrange(myset[m]['pos'], myset[m]['pos']+myset[m]['len']):
				try: mask[n] = 1
				except IndexError: pass

	for k in xrange(len(myset)):
		if myset[k]['chosen'] == 1:
			chosen.append(myset[k])
			#print '{},{},{},{},{},{},{},{},{}'.format(myset[k]['GQuadSeq'], myset[k]['g1Start'], myset[k]['g2Start'], myset[k]['g3Start'], myset[k]['g4Start'], myset[k]['gEnd'], myset[k]['glen'], myset[k]['len'], myset[k]['nscore'])

	#print 'Total: {}, Chosen: {}'.format(count, chosencnt)
	return chosen, count

def findall_QGRS(fa, output, MAX_GRS_LENGTH=45, GGROUPSIZE=2, LOOP_MIN=0, LOOP_MAX=36):
	
	uniq_cnt = 0
	overlap_cnt = 0
	search_space = 0
	body = ""

	for record in yield_fasta(fa):
		search_space = search_space+len(record.seq)
		chosenQGRS, nOverlap = QGRS(record.seq, MAX_GRS_LENGTH, GGROUPSIZE, LOOP_MIN, LOOP_MAX)
		uniq_cnt = uniq_cnt+len(chosenQGRS)
		overlap_cnt = overlap_cnt+nOverlap
		chosenQGRS = sorted(chosenQGRS, key=lambda k: k['g1Start'])
		for q in chosenQGRS:
			body += '{}\t{}\t{}\t{}\t{}\n'.format(record.ID, q['g1Start'], q['len'], q['nscore'], q['GQuadSeq'])

	summary = '#QGRS found: {}\n#QGRS found (including overlaps): {}\n#uniq/kb: {:.3f}\n#overlapping/kb: {:.3f}'.format(uniq_cnt, overlap_cnt, uniq_cnt*(1000.0/search_space), overlap_cnt*(1000.0/search_space))
	if output:
		with open(output, 'w') as fout:
			fout.write( '#Search parameters: QGRS Max Length: {}; Min G-Group Size: {}; Loop size: from {} to {}\n'.format(MAX_GRS_LENGTH, GGROUPSIZE, LOOP_MIN, LOOP_MAX) )
			fout.write( summary+"\n" )
			fout.write( '#ID\tPosition\tLength\tG-score\tQGRS\n' )
			fout.write( body )

	return uniq_cnt, overlap_cnt, uniq_cnt*(1000.0/search_space), overlap_cnt*(1000.0/search_space), search_space

if __name__ == '__main__':
	
	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--fasta', help="Nucleotide sequences in fasta format.", required=True)
	parser.add_argument('-o', '--output', help="QGRS output file.", required=False)
	parser.add_argument('--MAX-GRS-LENGTH', help="Maximum QGRS size. Default is 45.", default=45, required=False, type=int, choices=range(10, 46), metavar="[10-45]")
	parser.add_argument('--GGROUPSIZE', help="Minimum required G's. Default is 2", default=2, required=False, type=int, choices=range(2, 7), metavar="[2-6]")
	parser.add_argument('--LOOP_MIN', help="Minimum loop size. Default is 0.", default=0, required=False, type=int, choices=range(0, 37), metavar="[0-36]")
	parser.add_argument('--LOOP_MAX', help="Maximum loop size. Default is 36.", default=36, required=False, type=int, choices=range(0, 37), metavar="[0-36]")
	args = parser.parse_args()
	
	findall_QGRS(args.fasta, args.output, args.MAX_GRS_LENGTH, args.GGROUPSIZE, args.LOOP_MIN, args.LOOP_MAX)
