#!/usr/bin/env python

"""
Draw a VennDiagram based on a maximum of 5 transcript or gene name lists.
The VennDiagram is drawn by R.
"""

import argparse
import subprocess
from itertools import combinations
import os

def draw_venn_diagram(files, output_file="./venn_diagram.pdf", \
	labels=['A', 'B', 'C', 'D', 'E'], \
	colors=["#1492CC", "#F71E50", "#DFEDEF", "#00374C", "#350125"], \
	force_gene=False):

	N = len(files)
	r = range(N)

	if N > 1 and N < 6:

		# Create output_dir if it not exists
		output_dir = '/'.join(output_file.split('/')[:-1])
		if not os.path.exists(output_dir): os.makedirs(output_dir)

		# Parse transcripts/genes
		if force_gene: s = [ set([ line.rstrip().split('\t')[0].split('.')[0].upper() for line in open(f) if not line.startswith('#') ]) for f in files ]
		else: s = [ set([ line.rstrip().split('\t')[0].upper() for line in open(f) if not line.startswith('#') ]) for f in files ]

		with open("draw_venn.R", 'w') as fout:

			fout.write( "#!/usr/bin/env Rscript\n" )
			fout.write( "suppressMessages(library(VennDiagram))\n" )
			if N == 2: fout.write( "venn.plot <- draw.pairwise.venn(\n" )
			elif N == 3: fout.write( "venn.plot <- draw.triple.venn(\n" )
			elif N == 4: fout.write( "venn.plot <- draw.quad.venn(\n" )
			elif N == 5: fout.write( "venn.plot <- draw.quintuple.venn(\n" )

			for i in xrange(N): fout.write( "\tarea{} = {},\n".format(i+1, len(s[i])) )
			for i in xrange(2, N+1):
				for c in combinations(r, i):

					inter = s[c[0]].intersection(s[c[1]])
					if len(c) > 2:
						for j in xrange(2, len(c)):
							inter = inter.intersection(s[c[j]])

					if N == 2: fout.write( "\tcross.area = {},\n".format( len(inter) ) )
					else: fout.write( "\tn{} = {},\n".format( ''.join(str(n+1) for n in c), len(inter) ) )

			fout.write( "\tcategory = c({}),\n".format( ', '.join([ '"'+l+'"' for l in labels[:N] ]) ) )
			fout.write( "\tfill = c({}),\n".format( ', '.join([ '"'+c+'"' for c in colors[:N] ]) ) )
			fout.write( "\teuler.d = TRUE\n)\n" )
			fout.write( "pdf(file='{}{}');\n".format( output_file, '' if output_file.endswith('.pdf') else '.pdf' ) )
			fout.write( "grid.draw(venn.plot);" )

		subprocess.call("Rscript draw_venn.R", shell=True)
		subprocess.call("rm draw_venn.R Rplots.pdf", shell=True)

	else: print "ERROR: Too few or too many input files. Please supply 2 to 5 files."

if __name__ == '__main__':

	parser = argparse.ArgumentParser(description=__doc__)
	parser.add_argument('-f', '--files', required=True, nargs='+', help="Transcript/Gene lists")
	parser.add_argument('-o', '--output', default="./venn_diagram.pdf", help="Output path/file")
	parser.add_argument('-l', '--labels', nargs='+', default=['A', 'B', 'C', 'D', 'E'], help="Condition labels")
	parser.add_argument('-c', '--colors', nargs='+', default=["#1492CC", "#F71E50", "#DFEDEF", "#00374C", "#350125"], help="Condition colors")
	parser.add_argument('--FORCE-GENE', help="Force identifiers to gene format", action="store_true")
	args = parser.parse_args()

	draw_venn_diagram(args.files, args.output, args.labels, args.colors, args.FORCE_GENE)