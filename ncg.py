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
		self.centroid = self.getCentroid()

	def dist(self, a1, a2):
		return math.sqrt((a1.x-a2.x)*(a1.x-a2.x) + (a1.y-a2.y)*(a1.y-a2.y) + (a1.z-a2.z)*(a1.z-a2.z))

	# get current group geometric center
	def getCentroid(self):
		x=0.0
		y=0.0
		z=0.0
		n = len(self.atoms)
		for a in atoms:
			x+=atoms.x
			y+=atoms.y
			z+=atoms.z
		return (x/n, y/n, z/n)

	# grow contact group by comparing against input atoms
	# absorb the nearest one
	def grow(self, atoms):
		while len(self.atoms) < self.size:
			mindist = 999 # large enough
			for a in atoms:
				if a in self.atoms:
					continue
				if dist(a, self.centroid) < mindist:
					candidate = a
			self.atoms.append(a)

	def outStr(self):
		return '-'.joint(['%d%s' % (a.resSeq, a.chainID) for a in self.atoms])