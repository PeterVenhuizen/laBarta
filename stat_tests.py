#!/usr/bin/env python

from __future__ import division

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

def mc_perm_test(x, y, nmc=1000):
	"""
		Monte Carlo implementation of the permutation test
		
		http://biol09.biol.umontreal.ca/PLcourses/Statistical_tests.pdf
		http://stackoverflow.com/questions/24795535/pythons-implementation-of-permutation-test-with-permutation-number-as-input
		http://www.clayford.net/statistics/tag/permutation-tests/
	"""

	from itertools import combinations
	import numpy as np
	import math
	
	# Calculate the number combinations
	n = len(x) + len(y)
	r = len(x)
	N = math.factorial(n) / math.factorial(r) / math.factorial((n-r))

	# Do nmc iterations or all combinations when N < nmc
	k = 0
	diff = np.abs(np.mean(x) - np.mean(y))

	if N > nmc:

		# Do nmc random permutations
		for i in xrange(nmc):
			xy = np.concatenate([x, y])
			np.random.shuffle(xy)
			k += diff < np.abs(np.mean(xy[:n]) - np.mean(xy[n:]))
		return k / nmc

	else:

		# Check all permutations
		xy = np.concatenate([x, y])
		xy_indx = set([ i for i in xrange(n) ])
		for comb in combinations(xy_indx, r):
			new_x = [ xy[i] for i in comb ]
			new_y = [ xy[i] for i in set(xy_indx).difference(set(comb)) ]
			k += diff < np.abs(np.mean(new_x) - np.mean(new_y))
		return k / N

def fdr_correction(pvalues, rm_nan=True, fdr_method='fdr_bh'):
	'''
		Perform fdr correction with the possibility of ignoring NaN values.
		https://github.com/statsmodels/statsmodels/issues/2899
	'''

	from statsmodels.sandbox.stats.multicomp import multipletests
	import numpy as np

	if rm_nan:
		
		p = np.array(pvalues)
		mask = np.isfinite(p)
		pval_corrected = np.full(p.shape, np.nan)
		pval_corrected[mask] = multipletests(p[mask], method=fdr_method)[1]
		return pval_corrected

	else: 
		return multipletests(p[mask], method=fdr_method)[1]