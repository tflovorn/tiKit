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
#
# strategy = 'SquareSum' or 'Sum' - are PROCAR values squared before
# summation, or not?
def MarkSurfaceStates(procar, depth, threshold, strategy='SumSquare'):
    for kPoint in procar.kPoints:
        for band in kPoint.bands:
            marked = False
            for table in band.tables:
                isSurf, top, bottom = IsSurface(table, depth, threshold, strategy)
                band.top.append(top)
                band.bottom.append(bottom)
                if isSurf:
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
def IsSurface(table, depth, threshold, strategy):
    if strategy not in ['SumSquare', 'Sum']:
        print("error: invalid strategy")
        return False, 0.0, 0.0
    Ni = len(table.ions)
    sumTop, sumBottom, sumAll = 0.0, 0.0, 0.0
    # iterate over ions close to the top/bottom
    for i in range(1, Ni+1): 
        if strategy == 'SumSquare':
            sumAll += table.Ion(i).SquareSum()
        if i <= depth:
            # only here if i = 1, 2, ..., depth
            if strategy == 'SumSquare':
                sumTop += table.Ion(i).SquareSum()
                sumBottom += table.Ion(Ni - i + 1).SquareSum()
            elif strategy == 'Sum':
                sumTop += table.Ion(i).tot
                sumBottom += table.Ion(Ni - i + 1).tot
    if strategy == 'Sum':
        sumAll = abs(table.tot.tot)
    # surface weight above threshold?
    # note that sum(Top, Bottom, All) >= 0
    if sumAll < 1e-9:
        # total weight too small; assume not surface state
        return False, 0.0, 0.0
    if abs(sumTop/sumAll) > threshold or abs(sumBottom/sumAll) > threshold:
        # seems to be a surface state
        print('top: ' + str(sumTop) + ' bottom: ' + str(sumBottom) + ' total: ' + str(sumAll))
        return True, abs(sumTop/sumAll), abs(sumBottom/sumAll)
    # not above threshold
    return False, abs(sumTop/sumAll), abs(sumBottom/sumAll)

if __name__ == "__main__":
    # test - TODO arguments?
    with open('PROCAR', 'r') as f:
        procar = parseProcar.PROCAR(f, nonCol=True, lmDecomposed=True, storeIds=True)
        MarkSurfaceStates(procar, 3, 0.99)
