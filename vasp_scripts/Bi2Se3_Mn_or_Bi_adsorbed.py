#!/usr/bin/python
#
# Generate POSCAR for Bi2Se3-MnBi system with single Mn or Bi adsorbed.
# Hexagonal Bi2Se3 cells scaled in-plane as specified in the first argument.
#
# Usage: python Bi2Se3_MnBi_poscar.py (number of Bi/Se atoms per layer) (MnBi interface layer: 'Mn' or 'Bi') (MnBi interface position: 'aligned' or 'offset') (vacuum thickness in A)#
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
import sys, math

usageStr = "Usage: python Bi2Se3_MnBi_poscar.py (number of Bi/Se atoms per layer) (MnBi interface layer: 'Mn' or 'Bi') (MnBi interface position: 'aligned' or 'offset') (vacuum thickness in A)"

if len(sys.argv) != 5:
    print(usageStr)
    sys.exit(2)

# one Bi2Se3 quintuple-layer (QL) and one MnBi layer
N5L = 1
NMnBiL = 1

# number of atoms per layer in Bi2Se3 QL
scale = int(sys.argv[1])
if scale < 1:
    print(usageStr)
    print("Number of Bi/Se atoms per layer must be at least 1.")
    sys.exit(2)

# interface layer type: Bi-Se-Mn-Bi or Bi-Se-Bi-Mn
interface_layer_type = sys.argv[2]
if interface_layer_type not in ['Mn', 'Bi']:
    print(usageStr)
    print("Interface type must be Mn or Bi")
    sys.exit(2)

# interface layer position: aligned (A-A) or offset (A-B)
interface_layer_position = sys.argv[3]
if interface_layer_position not in ['aligned', 'offset']:
    print(usageStr)
    print("Interface type must be aligned (A-A) or offset (A-B)")
    sys.exit(2)

# vacuum thickness in Angstroms (absolute units)
c_vac_abs = float(sys.argv[4])
if c_vac_abs < 0.0:
    print(usageStr)
    print("Negative vacuum thickness not allowed.")
    sys.exit(2)

# lattice parameters in Angstroms (absolute units)
a_hex_abs = 4.138
c_1_abs = 1.7644
c_2_abs = 1.8799
c_vdw_abs = 2.2328
c_ql_abs = 2.0*c_1_abs + 2.0*c_2_abs
c_t_abs = c_ql_abs + c_vdw_abs
c_hex_abs = 3.0*c_t_abs

# MnBi lattice parameters in Angstroms (absolute units)
MnBi_a_hex_abs = 4.2827
MnBi_c_hex_abs = 6.1103

# assume interface distance = c_vdw
c_tot_abs = N5L*c_t_abs + 0.25*NMnBiL*MnBi_c_hex_abs + c_vac_abs

# c axis length in units relative to scale*a_hex_abs
c_tot = c_tot_abs / (scale*a_hex_abs)

# c axis lattice parameters in units relative to c_tot_abs
c_1 = c_1_abs / c_tot_abs
c_2 = c_2_abs / c_tot_abs
c_vdw = c_vdw_abs / c_tot_abs
c_ql = c_ql_abs / c_tot_abs
c_t = c_t_abs / c_tot_abs
c_hex = c_hex_abs / c_tot_abs
MnBi_c_hex = MnBi_c_hex_abs / c_tot_abs

# substitute in-plane distances relative to scale
subs_2_3 = 0.666667 / float(scale)
subs_1_3 = 0.333333 / float(scale)

def main():
    f = open('POSCAR', 'w')
    try:
        writeHexPOSCAR(f)
    finally:
        f.close()

def writeHexPOSCAR(f):
    header = "Bi2Se3/MnBi hex cell: 1x5 Bi2Se3 layers, in-plane scale " + str(scale) + "; 1 MnBi layers; Se-" + interface_layer_type + " interface; " + interface_layer_position + " interface position; vacuum " + str(c_vac_abs) + " Angstroms\n"
    f.write(header)
    f.write(str(float(scale)*a_hex_abs) + "\n")          # lattice parameter (scale)
    f.write("1.0 0.0 0.0\n")            # a_1
    f.write("-0.5 0.86602540378443864676 0.0\n")       # a_2
    f.write("0.0 0.0 " + str(c_tot) + "\n")    # a_3

    NMn, NBi = numMnBi()
    if NMn == 0:
        f.write(str(2*N5L*(scale**2) + NBi) + " " + str(3*N5L*(scale**2)) + "\n") # number of Bi/Se atoms
    else:
        f.write(str(2*N5L*(scale**2) + NBi) + " " + str(3*N5L*(scale**2)) + " " + str(NMn) + "\n") # number of Bi/Se/Mn atoms
    f.write("Direct\n") # direct coordinates (not cartesian)

    # atom positions: Bi in Bi2Se3 layers
    write_layer(f, "B", "Bi1", c_1)
    write_layer(f, "A", "Bi1p", c_1 + 2.0*c_2)

    # atom positions: Bi in MnBi layers
    if interface_layer_type == "Bi":
        if interface_layer_position == "aligned":
            f.write(str(subs_2_3) + " " + str(subs_1_3) + " " + str(c_t) + " Bi-B\n")
        else:
            f.write(str(subs_1_3) + " " + str(subs_2_3) + " " + str(c_t) + " Bi-C\n")

    # atom positions: Se in Bi2Se3 layers
    write_layer(f, "A", "Se1", 0.0)
    write_layer(f, "C", "Se2", c_1 + c_2)
    write_layer(f, "B", "Se1p", c_ql)

    # atom positions: Mn in MnBi layers
    if interface_layer_type == "Mn":
        if interface_layer_position == "aligned":
            f.write(str(subs_2_3) + " " + str(subs_1_3) + " " + str(c_t) + " Mn1-B\n")
        else:
            f.write(str(subs_1_3) + " " + str(subs_2_3) + " " + str(c_t) + " Mn1-C\n")

# Calculate the number of Mn and Bi atoms in MnBi.
def numMnBi():
    if interface_layer_type == "Mn":
        return 1, 0
    else:
        return 0, 1

# Create a layer of N_atoms = scale**2 with the given alignment
# ("A", "B", or "C") at the given c-axis position.
def write_layer(f, align, atom, c):
    subs_1 = 1.0 / float(scale)
    basis = []
    for i in range(scale):
        for j in range(scale):
            i_base = float(i)*subs_1
            j_base = float(j)*subs_1
            if align == "A":
                basis.append([i_base, j_base])
            elif align == "B":
                basis.append([i_base + subs_2_3, j_base + subs_1_3])
            else:
                basis.append([i_base + subs_1_3, j_base + subs_2_3])
    for pos in basis:
        f.write(str(pos[0]) + " " + str(pos[1]) + " " + str(c) + " " + atom + "-" + align + "\n")

if __name__ == "__main__":
    main()
