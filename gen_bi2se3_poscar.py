import sys

if len(sys.argv) != 3:
    print("usage: python gen_bi2se3_poscar.py N3QL c_vac\n")
    sys.exit(2)

# number of 15-layer cells (3x quintuple layers)
N3QL = int(sys.argv[1])
if N3QL < 1:
    print("usage: python gen_bi2se3_poscar.py N3QL c_vac\n")
    print("N3QL must be at least 1\n")
    sys.exit(2)

# vacuum length in Angstroms (absolute units)
c_vac_abs = float(sys.argv[2])

# lattice parameters in Angstroms (absolute units)
a_hex_abs = 4.138
c_1_abs = 1.7644
c_2_abs = 1.8799
c_vdw_abs = 2.2328
c_ql_abs = 2.0*c_1_abs + 2.0*c_2_abs
c_hex_abs = 3.0*(c_ql_abs + c_vdw_abs)
c_tot_abs = N3QL*c_hex_abs + c_vac_abs

# c axis length in units relative to a_hex_abs
c_tot = c_tot_abs / a_hex_abs

# c axis lattice parameters in units relative to c_tot_abs
c_1 = c_1_abs / c_tot_abs
c_2 = c_2_abs / c_tot_abs
c_vdw = c_vdw_abs / c_tot_abs
c_ql = c_ql_abs / c_tot_abs
c_hex = c_hex_abs / c_tot_abs

f = open('POSCAR', 'w')
try:
    header = "Bi2Se3: " + str(N3QL) + "x15 layers; vacuum " + str(c_vac_abs) + " Angstroms\n"
    f.write(header)
    f.write(str(a_hex_abs) + "\n")          # lattice parameter (scale)
    f.write("1.0 0.0 0.0\n")            # a_1
    f.write("-0.5 0.86602540378443864676 0.0\n")       # a_2
    f.write("0.0 0.0 " + str(c_tot) + "\n")    # a_3
    f.write(str(6*N3QL) + " " + str(9*N3QL) + "\n") # number of Bi/Se atoms
    f.write("Direct\n") # direct coordinates (not cartesian)
    # atom positions
    for n in range(N3QL): # A1, B1, C1
        f.write("0.666667 0.333333 " + str(n*c_hex + c_1) + " QL1-Bi1-B\n")
        f.write("0.0 0.0 " + str(n*c_hex + c_1 + 2.0*c_2) + " QL1-Bi1p-A\n")
        f.write("0.0 0.0 " + str(n*c_hex + c_vdw + c_ql + c_1) + " QL2-Bi1-A\n")
        f.write("0.333333 0.666667 " + str(n*c_hex + c_vdw + c_ql + c_1 + 2.0*c_2) + " QL2-Bi1p-C\n")
        f.write("0.333333 0.666667 " + str(n*c_hex + 2.0*c_vdw + 2.0*c_ql + c_1) + " QL3-Bi1-C\n")
        f.write("0.666667 0.333333 " + str(n*c_hex + 2.0*c_vdw + 2.0*c_ql + c_1 + 2.0*c_2) + " QL3-Bi1p-B\n")
        f.write("0.0 0.0 " + str(n*c_hex) + " QL1-Se1-A\n")
        f.write("0.333333 0.666667 " + str(n*c_hex + c_1 + c_2) + " QL1-Se2-C\n")
        f.write("0.666667 0.333333 " + str(n*c_hex + c_ql) + " QL1-Se1p-B\n")
        f.write("0.333333 0.666667 " + str(n*c_hex + c_vdw + c_ql) + " QL2-Se1-C\n")
        f.write("0.666667 0.333333 " + str(n*c_hex + c_vdw + c_ql + c_1 + c_2) + " QL2-Se2-B\n")
        f.write("0.0 0.0 " + str(n*c_hex + c_vdw + 2.0*c_ql) + " QL2-Se1p-A\n")
        f.write("0.666667 0.333333 " + str(n*c_hex + 2.0*c_vdw + 2.0*c_ql) + " QL3-Se1-B\n")
        f.write("0.0 0.0 " + str(n*c_hex + 2.0*c_vdw + 2.0*c_ql + c_1 + c_2) + " QL3-Se2-A\n")
        f.write("0.333333 0.666667 " + str(n*c_hex + 2.0*c_vdw + 3.0*c_ql) + " QL3-Se1p-C\n")
finally:
    f.close()
