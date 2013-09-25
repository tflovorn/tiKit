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
import sys
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import parseProcar, surface

def PlotBands(procar):
    # Iterate over bands to collect data.
    for i in range(1, procar.Nb+1):
        # Iterate over k-points to build (kx, energy) pairs for this band.
        # Entries in surface are 0 (bulk states) or 1 (surface states).
        # surface controls the color of the markers.
        kx, energy, surface = [], [], []
        for k in procar.kPoints:
            kx.append(k.kx)
            b = k.Band(i)
            energy.append(b.energy)
            if hasattr(b, 'surface') and b.surface:
                print(str(k.kx) + " " + str(i))
                surface.append(cm.gray(0.0))
            else:
                surface.append(cm.gray(0.99))
        # plot band as thin green line (to clarify connection btwn markers)
        plt.plot(kx, energy, 'g-', linewidth=1)
        # plot band as markers: bulk -> white, surface -> black
        plt.scatter(kx, energy, facecolors=surface, edgecolor='b',
                    marker='o', cmap=cm.gray, linewidth=1, s=50)

    # Display the plot
    plt.show()

if __name__ == "__main__":
    procarFileName = 'PROCAR'
    if len(sys.argv) > 1:
        procarFileName = sys.argv[1]
    procar = None
    with open(procarFileName, 'r') as f:
        procar = parseProcar.PROCAR(f, nonCol=True, lmDecomposed=True, storeIds=True)
    surface.MarkSurfaceStates(procar, 5, 0.70)
    PlotBands(procar)
