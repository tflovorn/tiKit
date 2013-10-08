#!/usr/bin/python
#
# Generate POSCAR for Bi2Se3-MnBi system. Hexagonal Bi2Se3 cells; possible
# to produce a cell smaller than the unit cell (which is 3 of the
# 5-layer cells).
# Usage: python Bi2Se3_MnBi_poscar.py (number of 5-layer Bi2Se3 cells) (number of MnBi layers) (MnBi interface layer: "Mn" or "Bi") (MnBi interface position: "A" or "B") (vacuum thickness in A)
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
import sys, math

usageStr = "Usage: python Bi2Se3_MnBi_poscar.py (number of 5-layer Bi2Se3 cells) (number of MnBi layers) (MnBi interface layer: 'Mn' or 'Bi') (MnBi interface position: 'aligned' or 'offset') (vacuum thickness in A)"

if len(sys.argv) != 6:
    print(usageStr)
    sys.exit(2)

# number of 5-layer Bi2Se3 cells (quintuple layers)
N5L = int(sys.argv[1])
if N5L < 1:
    print(usageStr)
    print("Number of 5-layer cells must be at least 1.")
    sys.exit(2)

# number of MnBi layers
NMnBiL = int(sys.argv[2])
if NMnBiL < 1:
    print(usageStr)
    print("Number of MnBi layers must be at least 1.")
    sys.exit(2)

# interface layer type: Bi-Se-Mn-Bi or Bi-Se-Bi-Mn
interface_layer_type = sys.argv[3]
if interface_layer_type not in ['Mn', 'Bi']:
    print(usageStr)
    print("Interface type must be Mn or Bi")
    sys.exit(2)

# interface layer position: aligned (A-A) or offset (A-B)
interface_layer_position = sys.argv[4]
if interface_layer_position not in ['aligned', 'offset']:
    print(usageStr)
    print("Interface type must be aligned (A-A) or offset (A-B)")
    sys.exit(2)

# vacuum thickness in Angstroms (absolute units)
c_vac_abs = float(sys.argv[5])
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

# c axis length in units relative to a_hex_abs
c_tot = c_tot_abs / a_hex_abs

# c axis lattice parameters in units relative to c_tot_abs
c_1 = c_1_abs / c_tot_abs
c_2 = c_2_abs / c_tot_abs
c_vdw = c_vdw_abs / c_tot_abs
c_ql = c_ql_abs / c_tot_abs
c_t = c_t_abs / c_tot_abs
c_hex = c_hex_abs / c_tot_abs
MnBi_c_hex = MnBi_c_hex_abs / c_tot_abs

def main():
    f = open('POSCAR', 'w')
    try:
        writeHexPOSCAR(f)
    finally:
        f.close()

def writeHexPOSCAR(f):
    header = "Bi2Se3/MnBi hex cell: " + str(N5L) + "x5 Bi2Se3 layers; " + str(NMnBiL) + " MnBi layers; Se-" + interface_layer_type + " interface; " + interface_layer_position + " interface position; vacuum " + str(c_vac_abs) + " Angstroms\n"
    f.write(header)
    f.write(str(a_hex_abs) + "\n")          # lattice parameter (scale)
    f.write("1.0 0.0 0.0\n")            # a_1
    f.write("-0.5 0.86602540378443864676 0.0\n")       # a_2
    f.write("0.0 0.0 " + str(c_tot) + "\n")    # a_3

    NMn, NBi = numMnBi()
    if NMn == 0:
        f.write(str(2*N5L + NBi) + " " + str(3*N5L) + "\n") # number of Bi/Se atoms
    else:
        f.write(str(2*N5L + NBi) + " " + str(3*N5L) + " " + str(NMn) + "\n") # number of Bi/Se/Mn atoms
    f.write("Direct\n") # direct coordinates (not cartesian)
    # atom positions: Bi in Bi2Se3 layers
    ql_state = 1
    for n in range(N5L):
        offset = math.floor(n / 3.0) * c_hex
        if ql_state == 1:
            f.write("0.666667 0.333333 " + str(offset + c_1) + " QL1-Bi1-B\n")
            f.write("0.0 0.0 " + str(offset + c_1 + 2.0*c_2) + " QL1-Bi1p-A\n")
        elif ql_state == 2:
            f.write("0.0 0.0 " + str(offset + c_vdw + c_ql + c_1) + " QL2-Bi1-A\n")
            f.write("0.333333 0.666667 " + str(offset + c_vdw + c_ql + c_1 + 2.0*c_2) + " QL2-Bi1p-C\n")
        elif ql_state == 3:
            f.write("0.333333 0.666667 " + str(offset + 2.0*c_vdw + 2.0*c_ql + c_1) + " QL3-Bi1-C\n")
            f.write("0.666667 0.333333 " + str(offset + 2.0*c_vdw + 2.0*c_ql + c_1 + 2.0*c_2) + " QL3-Bi1p-B\n")
        ql_state = (ql_state % 3) + 1
    # atom positions: Bi in MnBi layers
    Bi_offset = N5L * c_t
    if interface_layer_type == "Mn":
        Bi_offset += MnBi_c_hex / 4.0
    alignBi = interfaceLayerAlignment("Bi")
    for n in range(NBi):
        if n % 2 == 0:
            if alignBi == "B":
                f.write("0.666667 0.333333 " + str(Bi_offset + 0.5*n*MnBi_c_hex) + " Bi1-B\n")
            if alignBi == "C":
                f.write("0.333333 0.666667 " + str(Bi_offset + 0.5*n*MnBi_c_hex) + " Bi1-C\n")
            if alignBi == "A":
                f.write("0.0 0.0 " + str(Bi_offset + 0.5*n*MnBi_c_hex) + " Bi1-A\n")
        else:
            if alignBi == "B":
                f.write("0.333333 0.666667 " + str(Bi_offset + 0.5*n*MnBi_c_hex) + " Bi1-C\n")
            if alignBi == "C":
                f.write("0.0 0.0 " + str(Bi_offset + 0.5*n*MnBi_c_hex) + " Bi1-A\n")
            if alignBi == "A":
                f.write("0.666667 0.333333 " + str(Bi_offset + 0.5*n*MnBi_c_hex) + " Bi1-B\n")
    # atom positions: Se in Bi2Se3 layers
    ql_state = 1
    for n in range(N5L):
        offset = math.floor(n / 3.0) * c_hex
        if ql_state == 1:
            f.write("0.0 0.0 " + str(offset) + " QL1-Se1-A\n")
            f.write("0.333333 0.666667 " + str(offset + c_1 + c_2) + " QL1-Se2-C\n")
            f.write("0.666667 0.333333 " + str(offset + c_ql) + " QL1-Se1p-B\n")
        elif ql_state == 2:
            f.write("0.333333 0.666667 " + str(offset + c_vdw + c_ql) + " QL2-Se1-C\n")
            f.write("0.666667 0.333333 " + str(offset + c_vdw + c_ql + c_1 + c_2) + " QL2-Se2-B\n")
            f.write("0.0 0.0 " + str(offset + c_vdw + 2.0*c_ql) + " QL2-Se1p-A\n")
        elif ql_state == 3:
            f.write("0.666667 0.333333 " + str(offset + 2.0*c_vdw + 2.0*c_ql) + " QL3-Se1-B\n")
            f.write("0.0 0.0 " + str(offset + 2.0*c_vdw + 2.0*c_ql + c_1 + c_2) + " QL3-Se2-A\n")
            f.write("0.333333 0.666667 " + str(offset + 2.0*c_vdw + 3.0*c_ql) + " QL3-Se1p-C\n")
        ql_state = (ql_state % 3) + 1
    # atom positions: Mn in MnBi layers
    Mn_offset = N5L * c_t
    if interface_layer_type == "Bi":
        Mn_offset += MnBi_c_hex / 4.0
    alignMn = interfaceLayerAlignment("Mn")
    for n in range(NMn):
        if alignMn == "B":
            f.write("0.666667 0.333333 " + str(Mn_offset + 0.5*n*MnBi_c_hex) + " Mn1-B\n")
        if alignMn == "C":
            f.write("0.333333 0.666667 " + str(Mn_offset + 0.5*n*MnBi_c_hex) + " Mn1-C\n")
        if alignMn == "A":
            f.write("0.0 0.0 " + str(Mn_offset + 0.5*n*MnBi_c_hex) + " Mn1-A\n")

# calculate number of Mn and Bi atoms in MnBi
def numMnBi():
    if NMnBiL % 2 == 0:
        NMn = NMnBiL / 2
        NBi = NMnBiL / 2
        return NMn, NBi
    # NMnBiL is odd if we get here
    if interface_layer_type == "Mn":
        NMn = (NMnBiL + 1) / 2
        NBi = (NMnBiL - 1) / 2
        return NMn, NBi
    else:
        NMn = (NMnBiL - 1) / 2
        NBi = (NMnBiL + 1) / 2
        return NMn, NBi

# Get the alignment ("A", "B", or "C") of the first layer in MnBi matching
# the provided layerType,
def interfaceLayerAlignment(layerType):
    SeLayers = ["B", "A", "C"]
    lastSe = SeLayers[(N5L - 1) % 3]
    if interface_layer_position == "aligned":
        if layerType == interface_layer_type:
            # Se A -> (Mn or Bi) A
            return lastSe
        elif interface_layer_type == "Mn":
            # Bi second layer: Mn A -> Bi B
            return SeLayers[(N5L - 2) % 3]
        else:
            # Mn second layer: Bi A -> Mn C
            return SeLayers[N5L % 3]
    else:
        if layerType == interface_layer_type:
            # Se A -> (Mn or Bi) B
            return SeLayers[(N5L - 2) % 3]
        elif interface_layer_type == "Mn":
            # Bi second layer: Mn B -> Bi C
            return SeLayers[N5L % 3]
        else:
            # Mn second layer: Bi B -> Mn A
            return lastSe

if __name__ == "__main__":
    main()
