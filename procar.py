# Parses PROCAR file
#
# Copyright (c) 2013 Tim Lovorn (tflovorn@crimson.ua.edu)
# Released under the MIT License
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
#     The above copyright notice and this permission notice shall be included in
#     all copies or substantial portions of the Software.
#
#     THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
#     IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#     FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#     AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#     LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
#     OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
#     THE SOFTWARE.
#
import re, collections

# Represents the data stored in a PROCAR file.
# Contains properties nonCol, Nk, Nb, Ni, and kPoints.
# kPoints is a list of KPoint objects.
class PROCAR(object):
    # Create PROCAR object by reading from the file-like object procarFile.
    # nonCol = true if this is a non-collinear calculation (2 spins: 4 spinor
    # components); otherwise nonCol = false.
    def __init__(self, procarFile, nonCol):
        self.nonCol = nonCol

        # Start at the beginning of the file
        try:
            procarFile.seek(0)
        except:
            # if we can't seek, assume we are at the beginning
            print("warning: Couldn't seek to start of PROCAR file")

        # Get global data
        procarFile.readline()                # discard line 1
        globalLine = procarFile.readline()   # globals are on line 2
        # (number of k-points, bands, and ions)
        self.Nk, self.Nb, self.Ni = map(int, re.findall(r'\d+', globalLine))
        procarFile.readline()       # discard line 3 (empty)
        # Iterate over k-points.
        # The next line read should be the first line of the first k-point entry.
        self.kPoints = []
        for kId in range(1, self.Nk+1):
            # Extract k-point data.
            self.kPoints.append(KPoint(procarFile, self, kId))
            # Advance to the next k-point.
            procarFile.readline() # empty line
   
# Represents the data for one k-point.
# Contains properties kId, kx, ky, kz, weight, and bands.
# bands is a list of Band objects.
class KPoint(object):
    # Extract k-point data from procarFile.
    # The next line read should be the first line of the k-point entry.
    def __init__(self, procarFile, procar, kId):
        self.kId = kId
        kHead = procarFile.readline()
        # fixed positions since numbers overlap when there is a minus sign
        self.kx, self.ky, self.kz = map(float, [kHead[18:29], kHead[29:40],
                                                kHead[40:51]])
        # weight is the last thing in the line, isolated by spaces
        self.weight = float(kHead.rstrip().split(' ')[-1])
        procarFile.readline() # empty line
        # iterate over bands
        self.bands = []
        for bandId in range(1, procar.Nb+1):
            # Extract band data.
            self.bands.append(Band(procarFile, procar, bandId))
            # Advance to the next band.
            procarFile.readline() # empty line

# Represents the data for one band belonging to a specific k-point.
# Contains properties bandId, energy, occ, tables.
# tables is a list of IonTable objects.
class Band(object):
    # Extract band data from procarFile.
    # The next line read should be the first line of the band entry.
    def __init__(self, procarFile, procar, bandId):
        self.bandId = bandId
        bandHead = procarFile.readline().strip().split()
        # energy is the fifth group in the line, isolated by spaces
        self.energy = float(bandHead[4])
        # occ is the last group in the line, isolated by spaces
        self.occ = float(bandHead[-1])
        # get ion tables
        procarFile.readline() # empty line
        procarFile.readline()  # skip line containing "ion   s   py  pz"...etc
        self.tables = []
        numTables = 1
        if procar.nonCol:
            numTables = 4
        for tableId in range(1, numTables+1):
            # Extract table data.
            self.tables.append(IonTable(procarFile, procar, tableId))

# Represents a table of ionic data belonging to a (k-point, band) pair.
# Contains properties tableId, ions, total.
# ions is a list of Ion objects and total is an Ion object containing the
# total values summed over all ions.
class IonTable(object):
    # Extract ion table data from procarFile.
    # The next line read should be the first line of the table, containing ion 1
    def __init__(self, procarFile, procar, tableId):
        self.tableId = tableId
        # read data for each ion
        self.ions = []
        for ionId in range(1, procar.Ni+1):
            self.ions.append(Ion(procarFile, ionId))
        # read total
        self.ions.append(Ion(procarFile, 0))

# Represents data for one ion, belonging to a (k-point, band, ionTable).
# Contains properties ionId, s, py, pz, px, dxy, dyz, dz2, dxz, dx2, tot.
class Ion(object):
    def __init__(self, procarFile, ionId):
        self.ionId = ionId
        # each entry in the line is always separated by spaces
        l = procarFile.readline().strip().split()
        self.s, self.py, self.pz, self.px = map(float, [l[1], l[2], l[3], l[4]])
        self.dxy, self.dyz, self.dz2, self.dxz = map(float, [l[5], l[6], l[7], l[8]])
        self.dx2, self.tot = map(float, [l[9], l[10]])

if __name__ == "__main__":
    # test - TODO arguments?
    with open('PROCAR', 'r') as procarFile:
        procar = PROCAR(procarFile, True)
        print(procar.Nk, procar.Nb, procar.Ni)
        # having to index by one less than id's is annoying
        # TODO - is there a simple solution?
        print(procar.kPoints[0].bands[0].tables[0].ions[19].tot)
