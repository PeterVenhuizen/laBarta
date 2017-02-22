#!/usr/bin/env python

def hypergeometric_test(x, M, n, N, pval_cutoff=0.05):
	''' 
		hypergeom.sf(x, M, n, N)
		x => Number of successes in subgroup
		M => Population size
		n => Number of successes in population (M)
		N => Sample size
		hypergeom.sf tests for the random chance that you find the same or higher
		number of successes as in the sample. 
		
		http://homepages.ulb.ac.be/~dgonze/TEACHING/hypergeom.pdf
	'''

	from scipy.stats import hypergeom
	pval = hypergeom.sf(x, M, n, N)
	if pval < pval_cutoff:
		print 'x: {} M: {} n: {} N: {} p-value: {:.7g}'.format(x, M, n, N, pval)