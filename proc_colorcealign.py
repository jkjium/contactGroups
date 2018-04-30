import sys

import commp as cp
from protein import protein
from alignflat import palign

def main():
	if len(sys.argv)<3:
		print 'Usage: python proc_colorcealign.py cathpair.flat identifier'
		exit(1)

	flatfile = sys.argv[1]
	identifier = sys.argv[2]
	outfile = '%s.pml' % flatfile

	'''
	colormap = ['[0.4023, 0, 0.1211]',
				'[0.5508, 0.0469, 0.1445]',
				'[0.6953, 0.0938, 0.1680]',
				'[0.7656, 0.2344, 0.2344]',
				'[0.8359, 0.3750, 0.3008]',
				'[0.8945, 0.5117, 0.4023]',
				'[0.9531, 0.6445, 0.5078]',
				'[0.9727, 0.7500, 0.6445]',
				'[0.9883, 0.8555, 0.7773]',
				'[0.9766, 0.9102, 0.9102]',
				'[0.9648, 0.9648, 0.9648]',
				'[0.8906, 0.9297, 0.9531]',
				'[0.8164, 0.8945, 0.9375]',
				'[0.6953, 0.8320, 0.9023]',
				'[0.5703, 0.7695, 0.8672]',
				'[0.4180, 0.6719, 0.8164]',
				'[0.2617, 0.5742, 0.7617]',
				'[0.1953, 0.4883, 0.7148]',
				'[0.1289, 0.3984, 0.6719]',
				'[0.0742, 0.2930, 0.5234]',
				'[0.0195, 0.1875, 0.3789]']
	'''
	colormap = ['[0.5625,0.3828,0.6641]',
				'[0.9766,0.6484,0.7266]',
				'[0.5758,0.6969,0.6148]',
				'[0.4109,0.5828,0.9375]',
				'[0.6445,0.5156,0.4688]',
				'[0.5531,0.2094,0.2836]',
				'[0.9438,0.8734,0.5828]',
				'[0.5117,0.5000,0.6133]']

	fout = open(outfile, 'w')
	fout.write('bg_color white\n')
	#fout.write('set ribbon_smooth, 2\n')
	#fout.write('set ribbon_width, 6\n')
	#fout.write('cartoon tube\n')
	#fout.write('set sphere_scale, .5\n')
	for c in xrange(0, len(colormap)):
		fout.write('set_color kc%d, %s\n' % (c, colormap[c] ))

	count =0
	with open(flatfile) as fp:
		for line in fp:
			pml = []

			pa = palign(line.strip())

			#print '%s:\n%s\n%s\n' % (pa.name, pa.seqA, pa.seqB)
			pdbs = pa.pairnames()

			p1 = protein(pdbs[0]+'.aln.pdb')
			rmap1 = cp.posmap(pa.seqA.upper(), p1.seq.upper())
			if len(rmap1) == 0:
				continue
			pml.append('load %s.aln.pdb' % pdbs[0])
			pml.append('color gray90, %s.aln' % pdbs[0])

			p2 = protein(pdbs[1]+'.aln.pdb')
			rmap2 = cp.posmap(pa.seqB.upper(), p2.seq.upper())
			if len(rmap2) == 0:
				continue

			pml.append('load %s.aln.pdb' % pdbs[1])
			pml.append('color gray90, %s.aln' % pdbs[1])

			pml.append('cartoon automatic')
			pml.append('as cartoon')
			#pml.append('show sphere, name ca')
			if len(p1.ca)==len(p1.resDict) and len(p2.ca)==len(p2.resDict):
				r1 = p1.ca
				r2 = p2.ca
			else:
				r1 = p1.atomsbygmcenter()
				r2 = p2.atomsbygmcenter()

			posset = pa.alnposlist()
			for k in xrange(0, len(posset)):
				p1resi = []
				p2resi = []
				color = int(1.0 * k * len(colormap) / len(posset))
				#print color, repr(posset[k])
				#print repr(posset[k])
				for s in posset[k]:
					#print 'rmap2[%d]: %d' %  (s, rmap2[s]),
					p1resi.append(str(r1[rmap1[s]].resSeq))
					p2resi.append(str(r2[rmap2[s]].resSeq))

				pml.append('color kc%d, %s.aln and resi %s' % (color, pdbs[0], '+'.join(p1resi)))
				pml.append('color kc%d, %s.aln and resi %s' % (color, pdbs[1], '+'.join(p2resi)))

			pml.append('zoom')
			pml.append('save %s.%s.pse' % (pa.name, identifier))
			pml.append('save %s.%s.png' % (pa.name, identifier))

			pml.append('delete all')

			fout.write('%s\n' % '\n'.join(pml))
			count+=1
			print '%d alns processed.' % count
	print 'save to: %s'	% outfile
	fout.close()

if __name__ == '__main__':
	main()
