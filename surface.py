# surface.py: Uses data from PROCAR to identify surface states
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
import parseProcar

# For each ion table in procar, if the quantities in the 'total' column are
# concentrated close enough to the surface then mark the table as a surface
# state by adding the property table.surface=True; otherwise set
# table.surface=False. For a given (kpoint, band) pair, if any tables have
# the surface property set to true, set band.surface=True; otherwise set
# band.surface=False.
#
# Whether or not an ion table represents a surface state is controlled by
# the depth and threshold parameters. Starting with ion #1 and moving inward
# to ion #(depth), add the values in the total column together and call this
# 'sumTop'. On the other end, add the total column values for ion #(Ni)
# to ion #(Ni-depth+1) and call this 'sumBottom'.
# (TODO: better to make depth fractional?)
# The value table.tot.tot contains the sum over all ions of the total column;
# call this value 'sum'.
# If |sumBottom|/|sum| > threshold or |sumTop|/|sum| > threshold, then the
# table represents a surface state.
#
# TODO it would be nice to have an automatic procedure for choosing good
# values for depth and threshold.
def MarkSurfaceStates(procar, depth, threshold):
    for kPoint in procar.kPoints:
        for band in kPoint.bands:
            marked = False
            for table in band.tables:
                if IsSurface(table, depth, threshold):
                    table.surface = True
                    marked = True
                else:
                    table.surface = False
            if marked:
                print("marked: kPoint " + str(kPoint.kId) + " band " + str(band.bandId))
                band.surface = True

# Return true if the given table meets the requirements for a surface state
# as described in the documentation for MarkSurfaceStates. Return false
# otherwise.
def IsSurface(table, depth, threshold):
    Ni = len(table.ions)
    sumTop, sumBottom = 0, 0
    tot = abs(table.tot.tot)
    # bail on tot == 0.0; could in principle miss some surface states which
    # have weight concentrated near the surface but zero total
    if tot < 1e-9:
        return False
    # iterate over ions close to the top/bottom
    for i in range(1, depth+1): # i = 1, 2, ..., depth
        sumTop += table.Ion(i).tot
        sumBottom += table.Ion(Ni - i + 1).tot
    # surface weight above threshold?
    if abs(sumTop)/tot > threshold or abs(sumBottom)/tot > threshold:
        return True
    # not above threshold
    return False

if __name__ == "__main__":
    # test - TODO arguments?
    with open('PROCAR', 'r') as f:
        procar = parseProcar.PROCAR(f, nonCol=True, lmDecomposed=False, storeIds=True)
        MarkSurfaceStates(procar, 3, 0.99)
