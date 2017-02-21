#!/usr/bin/env python

"""
A collection of miscellaneous useful functions.
"""

def roundrobin(*iterables):
    "roundrobin('ABC', 'D', 'EF') --> A D E B F C"
    # Recipe credited to George Sakkis
    pending = len(iterables)
    nexts = cycle(iter(it).next for it in iterables)
    while pending:
        try:
            for next in nexts:
                yield next()
        except StopIteration:
            pending -= 1
            nexts = cycle(islice(nexts, pending))

def get_levenshtein_distance(s1, s2):
	''' Determine the levenshtein distance of two sequence, codons, etc. '''
    
	if len(s1) > len(s2):
		s1, s2 = s2, s1
	distances = range(len(s1) + 1)
	for index2, char2 in enumerate(s2):
		newDistances = [index2+1]
		for index1, char1 in enumerate(s1):
			if char1 == char2:
				newDistances.append(distances[index1])
			else:
				newDistances.append(1 + min((distances[index1], distances[index1+1], newDistances[-1])))
		
		distances = newDistances
		
	return distances[-1]

def find_potential_PTC(seq, RF="+1"):
	''' Look for codons which could mutate into a stop codon by a SINGLE mutation. 
		Return a dictionary containing the codon positions and the codon itself. '''
	
	seq = seq.upper().replace('U', 'T')
	seq = get_frame_seq(seq, RF)
	seq_len = len(seq)
	stop_codons = ["TAG", "TAA", "TGA"]
	
	dangerous_codons = {}
	for i in xrange(0, seq_len, 3):
		codon = seq[i:i+3]
		for stop in stop_codons:
			if get_levenshtein_distance(codon, stop) == 1: dangerous_codons["%d-%d" % (i, i+3)] = codon
	
	return dangerous_codons

def iupac2regex(iupac_seq, nMatch=False):
	
	code = { 'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
		'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]', 'K': '[GT]', 'M': '[AC]',
		'B': '[CGT]', 'D': '[AGT]', 'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]' }
	
	if nMatch:
	
		reg = []
		score = 1
		for b in list(iupac_seq):
			m = code[b.upper()]
			reg.append(m)
			if '[' in m: score = score * (len(m)-2)
		
		return ''.join(reg), score

	else:
		return ''.join([ code[b.upper()] for b in list(iupac_seq) ])

def enumerate_iupac(iupac_seq, filename=''):
	
	IUPAC_translation = { 'A': 'A',
		'B': 'CGT',
		'C': 'C',
		'D': 'AGT',
		'G': 'G',
		'H': 'ACT',
		'I': 'ACGT',
		'K': 'GT',
		'M': 'AC',
		'N': 'ACGT',
		'R': 'AG',
		'S': 'CG',
		'T': 'T',
		'V': 'ACG',
		'W': 'AT',
		'Y': 'CT'
	}
	
	seqArray = []
	iupac_seq = iupac_seq.upper()
	for base in iupac_seq:
		bases = IUPAC_translation[base]
		nBases = len(bases)
		if nBases > 1:
			if len(seqArray) == 0:
				for b in bases: seqArray += [b]
			else:
				tmpArray = []
				for s in seqArray:
					for b in bases: tmpArray += [s+b]
				seqArray = tmpArray
		elif len(seqArray) > 0:
			tmpArray = []
			for s in seqArray: tmpArray += [s+bases]
			seqArray = tmpArray
		else:
			seqArray = [''.join(seqArray) + bases]
			
	if len(filename):
		with open(filename, 'w') as fout:
			c = 1
			for s in seqArray: 
				fout.write('>%d\n%s\n' % (c, s))
				c += 1
	else: 
		for s in seqArray: print s