import sys
import math
from atom import atom
class ncg(object):
	#init starts from a single atom
	#the final result(overlapped ncgs) rely on each atoms
	def __init__(self, atom, size):
		self.size = size
		self.atoms = []
		self.atoms.append(atom)
		self.centroid = (atom.x, atom.y, atom.z)

	def dist(self, a):
		return math.sqrt((a.x-self.centroid[0])*(a.x-self.centroid[0]) + (a.y-self.centroid[1])*(a.y-self.centroid[1]) + (a.z-self.centroid[2])*(a.z-self.centroid[2]))

	# get current group geometric center
	def calcCentroid(self):
		x=0.0
		y=0.0
		z=0.0
		n = len(self.atoms)
		for a in self.atoms:
			x+=a.x
			y+=a.y
			z+=a.z
		return (x/n, y/n, z/n)

	# grow contact group by comparing against input atoms
	# absorb the nearest one
	def grow(self, atoms):
		while len(self.atoms) < self.size:
			#print 'ncg::grow(): extracting %d order contact.' % len(self.atoms)
			mindist = 999 # large enough
			candidate = -1
			for a in atoms:
				if a in self.atoms:
					continue
				d = self.dist(a)
				if d < mindist:
					mindist = d
					candidate = a
			self.atoms.append(candidate)
			self.centroid = self.calcCentroid()

	def outStr(self):
		return ' '.join(['%s%d' % (a.chainID, a.resSeq) for a in self.atoms]) 