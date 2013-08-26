# Generate POSCAR for MnBi system
# Usage: python MnBiposcar.py (number of 4-layer cells) (vacuum thickness in A)
#
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
import sys

usageStr = "Usage: python MnBiposcar.py (number of 4-layer cells) (vacuum thickness in A)"

if len(sys.argv) != 3:
    print(usageStr)
    sys.exit(2)

# number of 4-layer cells
N4L = int(sys.argv[1])
if N4L < 1:
    print(usageStr)
    print("Number of 4-layer cells must be at least 1.")
    sys.exit(2)

# vacuum thickness in Angstroms (absolute units)
c_vac_abs = float(sys.argv[2])
if c_vac_abs < 0.0:
    print(usageStr)
    print("Negative vacuum thickness not allowed.")
    sys.exit(2)

# lattice parameters in Angstroms (absolute units)
a_hex_abs = 4.2827
c_hex_abs = 6.1103
c_tot_abs = N4L*c_hex_abs + c_vac_abs

# c axis length in units relative to a_hex_abs
c_tot = c_tot_abs / a_hex_abs

# c axis lattice parameters in units relative to c_tot_abs
c_hex = c_hex_abs / c_tot_abs

def main():
    f = open('POSCAR', 'w')
    try:
        writeHexPOSCAR(f)
    finally:
        f.close()

def writeHexPOSCAR(f):
    header = "MnBi hex cell: " + str(N4L) + "x4 layers; vacuum " + str(c_vac_abs) + " Angstroms\n"
    f.write(header)
    f.write(str(a_hex_abs) + "\n")          # lattice parameter (scale)
    f.write("1.0 0.0 0.0\n")            # a_1
    f.write("-0.5 0.86602540378443864676 0.0\n")       # a_2
    f.write("0.0 0.0 " + str(c_tot) + "\n")    # a_3
    f.write(str(2*N4L) + " " + str(2*N4L) + "\n") # number of Mn/Bi atoms
    f.write("Direct\n") # direct coordinates (not cartesian)
    # atom positions
    for n in range(N4L):
        f.write("0.0 0.0 " + str(float(n) * c_hex) + " Mn1-A\n")
        f.write("0.0 0.0 " + str((float(n)+0.5) * c_hex) + " Mn2-A\n")
    for n in range(N4L): 
        f.write("0.666667 0.333333 " + str((float(n)+0.25) * c_hex) + " Bi1-B\n")
        f.write("0.333333 0.666667 " + str((float(n)+0.75) * c_hex) + " Bi2-C\n")

if __name__ == "__main__":
    main()
