#!/usr/bin/python
## eigenval2foo
## Copyright (C) 2005 Crutcher Dunnavant <crutcher@gmail.com>
#            (C) 2013 Tim Lovorn <tflovorn@crimson.ua.edu>

## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

import os
import sys
import math
import string

# ====================================================================== #
# Status

def status(str):
	sys.stderr.write(str)
	sys.stderr.flush()

# ====================================================================== #
## EIGENVAL files have the following headers, which I only partially
## understand.
## UPDATED descriptions pulled from thread at:
## http://cms.mpi.univie.ac.at/vasp-forum/forum_viewtopic.php?4.7 ; post on Tue Oct 24 2006, 01:23PM
## ORIGINAL description: {?}  {?}  {?}  {num_spins:=(1|2)}
## UPDATED description: {number of atoms} {number of atoms (again)} {NBLOCK*KBLOCK} {1 = non-spin-polarized; 2 = spin-polarized}
#    2    2    1   1
## UPDATED description: {cell volume} {lattice vector a in m} {lattice vector b in m} {lattice vector c in m} {POTIM * 1e-15}
#  0.8323070E+01  0.2866000E-09  0.2866000E-09  0.2866000E-09  0.5000000E-15
## UPDATED description: {temperature} <-- TODO check this, EDIFF seems more likely
#   1.00000000000000E-004
#  CAR
# MgO 
##  ORIGINAL description: {?}  {num_k}  {num_bands}
##  UPDATED description: {number of electrons} {number of k points} {number of bands}
#    8   20    8

## Each k-point record looks like this:
# << Blank Line >>
## {a}            {b}            {c}            {weight}
#  0.0000000E+00  0.0000000E+00  0.0000000E+00  0.5000000E-01
## band: the band id for this data entry at this k-point
## point1: the data entry for the first spin.
## point2: the data entry for the second spin, if there is one.
## {band}  {point1}     [{point2}]
#   1      -11.9377
#   2        5.9077
#   3        5.9078
#   4        5.9078
#   5       11.8920
#   6       22.4112
#   7       22.4112
#   8       22.4112
 
class EIGENVAL:
	def __init__(self, path = None):
		self.load_from_path(path)

	def load_from_path(self, path):
		lines = open(path).readlines()

		hdr1 = lines[0].split()
		# Is this a 1 or 2 spin file?
		num_spins = int(hdr1[3])
		(junk, num_k, num_bands) = map(int, lines[5].split())

		def listoflists(n):
			return map(lambda x: [], range(n))

		if num_spins == 2:
			self.points = \
				[listoflists(num_bands),
				 listoflists(num_bands)]
		else:
			self.points = [listoflists(num_bands)]

		self.kpoints = []
		self.weights = []

		k_record_size = 2 + num_bands

		for k in range(num_k):
			# k-record structure:
			# 1 blank line.
			# 1 k-point line.
			# num_bands band lines.
			o = 6 + k * k_record_size # o = line number for start of current k-record

			a = map(float, lines[o + 1].split()) # o+1 --> k-point line (skip blank line)

			self.kpoints.append(tuple(a[:3]))
			self.weights.append(a[3])

			for b in range(num_bands):
				p = map(float, lines[o + 2 + b].split()[1:])

				for s in range(num_spins):
					self.points[s][b].append(p[s])

	def getNumSpins(self):
		return len(self.points)

	def getNumBands(self):
		return len(self.points[0])

	def getNumKpoints(self):
		return len(self.points[0][0])

	def projectKpoints(self):
		l = []
		for (a, b, c) in self.kpoints:
            # k = (kx, ky, kz) --> |k|
			l.append(math.sqrt(a**2 + b**2 + c**2)) 
		return l

	def __str__(self):
		info = {
			'num_spins' : self.getNumSpins(),
			'num_bands' : self.getNumBands(),
			'num_kpoints' : self.getNumKpoints(),
		}
		rstr = "    ?    ?    ?    %(num_spins)d\n" % info
		rstr += "  ??  ??  ??  ??  ??\n"
		rstr += "   ?\n"
		rstr += "  ???\n"
		rstr += " ???\n"
		rstr += "    ?    %(num_kpoints)d    %(num_bands)d\n" % info

		for k in range(self.getNumKpoints()):
			rstr += "\n  %E %E %E %E\n" % \
				(self.kpoints[k][0],
				 self.kpoints[k][1],
				 self.kpoints[k][2],
				 self.weights[k])

			for b in range(self.getNumBands()):
				rstr += "   %d" % (b + 1)
				for s in range(self.getNumSpins()):
					rstr += "      %f" % \
						self.points[s][b][k]
				rstr += "\n"
		return rstr

# ====================================================================== #
# Relax Curves

def D(L):
	# Given: [p1, p2, p3, ..., pn]
	# Return: [(p2 - p1), (p3 - p2), ..., (pn - pn-1)]
	return map(lambda p1, p2: p2 - p1, L[:-1], L[1:])
	
def cset_cost(cset):
	# {Assumption: delta-x is constant.}
	# Return the sum of the squares of the second-order deltas
	# of the curves in the cset.
	cost = 0
	for curve in cset:
		curve_cost = 0
		for d2 in D(D(curve)):
			curve_cost = curve_cost + d2**2
		cost = cost + curve_cost

	return cost

def permute_cset(cset, i, j):
	# Extract two cset.
	ab = cset[j]
	cd = cset[j + 1]

	# Cross the extracted cset
	ad = ab[:i] + cd[i:]
	cb = cd[:i] + ab[i:]

	# Return the new cset set.
	return cset[:j] + [ad, cb] + cset[j+2:]

def relax_cset(cset):
	# {Assumption: delta-x is constant.}

	cost = cset_cost(cset)

	M = len(cset)
	N = len(cset[0])

	c = 0
	updated = 1
	while updated:
		updated = 0
		c = c + 1

		status("=round %d: " % c)

		for j in range(M - 1):
			for i in range(N - 1):

				new_cset = permute_cset(cset, i, j)
				new_cost = cset_cost(new_cset)

				if new_cost < cost:
					cset = new_cset
					cost = new_cost
					updated = 1

					status("+")
				else:
					status(".")

		status("\n")

	status("=done.\n")

	return cset

# ====================================================================== #

usage = """usage: eigenval2foo [OPTIONS] INFILE FORMAT OUTFILE 

Performs various conversions on EIGENVAL files.

Available options:
	-r	Relax the curves.
	-1	Select only the first spin (default).
	-2	Select only the second spin.
	-0	Select all spins.
	(-0 is the same as -1 if there is only one spin in the input file).

Available formats:
	gnuplot
	csv

Examples:
	eigenval2foo -2 EIGENVAL csv EIGENVAL
		- Extracts the second spin from the file "EIGENVAL", and
		  writes it to the CSV data file "EIGENVAL.csv"

	eigenval2foo -r0 EIGENVAL gnuplot R
		- Extracts all spins from the file "EIGENVAL", relaxes
		  them, and writes the results to the gnuplot data file,
		  "R.dat". Generates the gnuplot script "R.plt", and runs
		  gnuplot on it, generating "R.png".
"""

def main():
	if len(sys.argv) == 1:
		print usage
		return

	# Defaults.
	relax = 0
	spin = 1

	for i in range(1, len(sys.argv)):
		arg = sys.argv[i]
		if arg[0] == '-':
			for c in arg[1:]:
				if c == "r":
					relax = 1
				elif c == "1":
					spin = 1
				elif c == "2":
					spin = 2
				elif c == "0":
					spin = 0
		else:
			break

	try:
		(in_path, format, out_name) = sys.argv[i:]
	except:
		print usage
		return

	status('Reading "%s"\n' % in_path)
	e = EIGENVAL(in_path)

	# Handle spin and relaxation interactions.
	if spin == 0 and e.getNumSpins() == 1:
		spin = 1

	if spin == 2 and e.getNumSpins() == 1:
		print 'Error, "%s" contains only one spin.' % in_path
		return

	if relax:
		if spin == 1:
			status('Relaxing Spin 1\n')
			A = relax_cset(e.points[0])

		elif spin == 2:
			status('Relaxing Spin 2\n')
			A = relax_cset(e.points[1])

		else:
			status('Relaxing Spin 1\n')
			A = relax_cset(e.points[0])
			
			status('Relaxing Spin 2\n')
			A = A + relax_cset(e.points[1])
	else:
		if spin == 1:
			A = e.points[0]

		elif spin == 2:
			A = e.points[1]

		else:
			A = e.points[0] + e.points[1]

	M = len(A)
	N = len(A[0])

#	X = e.projectKpoints()
	step = 1.0 / float(N - 1)
	X = map(lambda x, s=step: x * s, range(N))

	# A is a list of curves, which are lists of points.
    # ---> A follows the order of points in KPOINTS/EIGENVAL
	# X is a list of points for the X axis.
	# N = len(X) = len(A[i]), for all i
	# M = len(A)

	# Handle formats.
	if format == "gnuplot":
        # to clarify the zip operation:
        # >>> X = [0, 1, 2]
        # >>> [X]
        # [[0, 1, 2]]
        # >>> A = [[3, 4, 5], [6, 7, 8], [9, 10, 12]]
        # >>> apply(zip, [X] + A)
        # [(0, 3, 6, 9), (1, 4, 7, 10), (2, 5, 8, 12)]
        # >>> [X] + A
        # [[0, 1, 2], [3, 4, 5], [6, 7, 8], [9, 10, 12]]
		B = apply(zip, [X] + A)

		dat_path = out_name + '.dat'

		outf = open(dat_path, "w")
		status('Writting "%s"\n' % dat_path)
		for row in B:
			print >> outf, "\t".join(map(str, row))
		outf.close()

		plot_path = out_name + '.plt'
		status('Writting "%s"\n' % plot_path)
		poutf = open(plot_path, 'w')

	#	print >> poutf, "set terminal png"
	#	print >> poutf, "set output '%s.png'" % out_name
		out_pic_name = out_name + '.ps'
		print >> poutf, "set terminal postscript"
		print >> poutf, "set output '%s'" % out_pic_name

	#	print >> poutf, gnuplot_script_old
		print >> poutf, gnuplot_script

		print >> poutf, "unset key"

		poutf.write('plot ')
		l = map \
			(lambda x, d = dat_path:
				'"%s"u 1:%d title "%d"' % (d, x, x),
			 range(2, len(A) + 2))
		poutf.write(", ".join(l))
		poutf.write("\n")

		poutf.close()

		status('Generating "%s"\n' % out_pic_name)
		os.system('gnuplot "%s"' % plot_path)
		status('Generating "%s.pdf"\n' % out_name)
		os.system('ps2pdf "%s"' % out_pic_name)
		return

	elif format == "csv":
		B = apply(zip, [X] + A)

		import csv

		outf = open(out_name + '.csv', "w")
		status('Writting "%s"\n' % out_name)

		writer = csv.writer(outf)
		for row in B:
			writer.writerow(row)

		outf.close()
		return
	
	else:
		print "Error: Unknown format '%s'" % format
		return

# ====================================================================== #

gnuplot_script = r"""
set style data linespoints
set style function lines
"""
# ====================================================================== #

if __name__ == "__main__":
	main()

