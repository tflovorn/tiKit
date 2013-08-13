import sys, json, math
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt

# 3D model notation for Gamma matrices
Ident = np.array([[1, 0, 0, 0],
                  [0, 1, 0, 0],
                  [0, 0, 1, 0],
                  [0, 0, 0, 1]])
Gamma1 = np.array([[0, 0, 0, 1],
                   [0, 0, 1, 0],
                   [0, 1, 0, 0],
                   [1, 0, 0, 0]])
Gamma2 = np.array([[0, 0, 0, -1j],
                   [0, 0, -1j, 0],
                   [0, 1j, 0, 0],
                   [1j, 0, 0, 0]])
Gamma3 = np.array([[0, 1, 0, 0],
                   [1, 0, 0, 0],
                   [0, 0, 0, -1],
                   [0, 0, -1, 0]])
Gamma4 = np.array([[0, -1j, 0, 0],
                   [1j, 0, 0, 0],
                   [0, 0, 0, -1j],
                   [0, 0, 1j, 0]])
Gamma5 = np.array([[1, 0, 0, 0],
                   [0, -1, 0, 0],
                   [0, 0, 1, 0],
                   [0, 0, 0, -1]])

# Construct a Hermitian matrix contaning the upper triangular part given
# by M. The diagonal elements of M must be real and the elements below
# the diagonal must contain all zeros.
def makeHermitian(M):
    return M + M.T.conj() - np.diag(M.diagonal())

# Parse command-line arguments calcType, kpointsFileName, outFileName
def parseArgs():
    if len(sys.argv) != 4:
        print("Usage: ti3D_eigen calcType kpointsFileName outFileName")
        print("calcType should be 8band, 4band, or mnk12")
        sys.exit(2)
    return sys.argv[1], sys.argv[2], sys.argv[3]

# Parse kpoints file
def getKpoints(kpointsFileName):
    # KPOINTS format: ignore lines 0, 2, 3; line 1 = # of points in each range;
    # line 4-5: range 1; line 7-8: range 2; ...
    # k values returned are in 2pi/(lattice constant) units.
    kpointsFile = open(kpointsFileName, 'r')
    lines = kpointsFile.readlines()
    kpointsFile.close()

    kNum = int(lines[1])
    lineNum = 4
    kBounds = []
    while lineNum < len(lines):
        # take lineNum and lineNum+1 and convert each to a list of floats
        start = map(float, lines[lineNum].split(" "))
        stop = map(float, lines[lineNum+1].split(" "))
        kBounds.append([start, stop])
        lineNum += 3

    kpoints = []
    # iterate over (kStart, kStop) pairs and generate points
    for pair in kBounds:
        start, stop = np.array(pair[0]), np.array(pair[1])
        step = (stop - start)/(kNum - 1)
        # make kNum points between start and stop, including boundaries
        for ptIndex in range(kNum):
            point = start + ptIndex*step
            kpoints.append(point)
    return kpoints

# Return a Hamiltonian function (k -> H_k) of the type specified
def HamiltonianFn(calcType):
    fns = {"8band": Hamiltonian_8band, "4band": Hamiltonian_4band, "mnk12": Hamiltonian_mnk12}
    if calcType not in fns:
        print("Usage: ti3D_eigen calcType kpointsFileName outFileName")
        print("calcType should be 8band, 4band, or mnk12")
        sys.exit(2)
    with open(calcType + ".json", 'r') as propsFile:
        props = json.load(propsFile)
        return fns[calcType](props)

# Generate 8-band Hamiltonian function with properties given by p.
# This Hamiltonian is from Liu et al PRB 82, 045122 (2010).
# Complex numbers in p are encoded as "varName" + ("re" or "im").
def Hamiltonian_8band(p):
    k_plus = lambda k: k[0] + 1j*k[1]
    k_minus = lambda k: k[0] - 1j*k[1]
    f = lambda k, F, K: F * k[2]**2 + K * (k[0]**2 + k[1]**2)
    g = lambda k, U, V: U * k[2] * k_plus(k) + V * k_minus(k)**2

    def H(in_k):
        # Convert k from 2pi/(lattice vector) units to 1/A.
        # TODO: ensure that this assignment of a and c is correct.
        a_hex = 4.138
        c_hex = 28.64
        k = [0, 0, 0]
        k[0] = 2.0 * math.pi * in_k[0] / a_hex
        k[1] = 2.0 * math.pi * in_k[1] / a_hex
        k[2] = 2.0 * math.pi * in_k[2] / c_hex

        # Assume factors of (hbar)^2/(2m) and 2/hbar are absorbed into
        # constants. TODO: check that this is correct.
        diagonal = [f(k, p["F1"], p["K1"]), f(k, p["F1"], p["K1"]),
                    f(k, p["F3"], p["K3"]), f(k, p["F3"], p["K3"]),
                    f(k, p["F5"], p["K5"]), f(k, p["F5"], p["K5"]),
                    f(k, p["F7"], p["K7"]), f(k, p["F7"], p["K7"])]

        top = np.zeros([8, 8], dtype=np.complex128)
        top[0] = [0, 0, k[2]*p["Q1"], p["P1"]*k_minus(k), p["Q2"]*k_plus(k),
                  k_plus(k)*p["P2"], k[2]*p["Q3"], k_minus(k)*p["P3"]]
        # P1,2,3 and Q1,2,3 are real so ignore complex conjugate
        top[1] = [0, 0, k_plus(k)*p["P1"], -k[2]*p["Q1"], -p["P2"]*k_minus(k),
                  p["Q2"]*k_minus(k), p["P3"]*k_plus(k), -p["Q3"]*k[2]]

        U35 = p["U35re"] + 1j*p["U35im"]
        V35 = p["V35re"] + 1j*p["V35im"]
        U36, V36 = U35.conjugate(), -V35.conjugate()
        top[2] = [0, 0, 0, 0, g(k, U35, V35), g(k, U36, V36),
                  f(k, p["F37"], p["K37"]), 
                  -g(-k, p["U47"], 1j*p["V47im"]).conjugate()]
        top[3] = [0, 0, 0, 0, g(-k, U36, V36).conjugate(),
                  -g(-k, U35, V35).conjugate(), g(k, p["U47"], 1j*p["V47im"]),
                  f(-k, p["F37"], p["K37"]).conjugate()]

        U58 = p["U58re"] + 1j*p["U58im"]
        V58 = p["V58re"] + 1j*p["V58im"]
        U68, V68 = -U58.conjugate(), V58.conjugate()
        top[4] = [0, 0, 0, 0, 0, 0, -g(-k, U68, V68).conjugate(),
                  g(k, U58, V58)]
        top[5] = [0, 0, 0, 0, 0, 0, g(-k, U58, V58).conjugate(),
                  g(k, U68, V68)]

        return makeHermitian(top + np.diag(diagonal))

    return H

# Generate 4-band Hamiltonian function with properties given by p.
# This Hamiltonian is also from Liu 2010.
def Hamiltonian_4band(p):
    # (k_parallel)^2
    kp2 = lambda k: k[0]**2 + k[1]**2
    # k_z^2
    kz2 = lambda k: k[2]**2
    # functions in H0, Eq. 16, Liu 2010
    epsilon = lambda k: p["C0"] + p["C1"]*kz2(k) + p["C2"]*kp2(k)
    M = lambda k: p["M0"] + p["M1"]*kz2(k) + p["M2"]*kp2(k)
    
    A2 = -0.5 * p["A0"]
    if "A2" in p:
        A2 = p["A2"]
    A = lambda k: p["A0"] + A2*kp2(k)

    B2 = -0.5 * p["B0"]
    if "B2" in p:
        B2 = p["B2"]
    B = lambda k: p["B0"] + B2*kz2(k)

    # parenthetic expression in H3, Eq. 17, Liu 2010
    q = lambda kx, ky: kx**3 - 3.0*kx*(ky**2)

    def H(in_k):
        # Convert k from 2pi/(lattice vector) units to 1/A.
        # TODO: ensure that this assignment of a and c is correct.
        a_hex = 4.138
        c_hex = 28.64
        k = [0, 0, 0]
        k[0] = 2.0 * math.pi * in_k[0] / a_hex
        k[1] = 2.0 * math.pi * in_k[1] / a_hex
        k[2] = 2.0 * math.pi * in_k[2] / c_hex

        H0 = epsilon(k)*Ident + M(k)*Gamma5 + B(k)*Gamma4*k[2] + A(k)*(Gamma1*k[1] - Gamma2*k[0])
        H3 = p["R1"]*Gamma3*q(k[0], k[1]) - p["R2"]*Gamma4*q(k[1], k[0])
        return H0 + H3

    return H

# Generate mnk12 Hamiltonian function with properties given by p
def Hamiltonian_mnk12(p):
    d = lambda k: p["M"] - 2.0*p["B"] + 2.0*p["B"] * (math.cos(k[0])
            + math.cos(k[1]) - 2.0)
    def H(k):
        # Convert k from 2pi/(lattice vector) units to 1/(lattice vector).
        # Lattice vector factor eliminated in sin(kx * a), etc.
        k = 2.0 * math.pi * k

        diagonal = np.diag([p["C"]]*4)
        sin_part = p["A"]*(Gamma2 * math.sin(k[0]) + Gamma1 * math.sin(k[1])
                           + Gamma4 * math.sin(k[2]))
        cos_part = Gamma5 * (2.0*p["B"] * math.cos(k[2]) + d(k))
        return diagonal + sin_part + cos_part

    return H

# Write the output for one kpoint
def writeOutput(k, eigenvals, eigenkets, outFile):
    #TODO
    outFile.write(str(k) + "\n")
    outFile.write(str(eigenvals) + "\n")
    outFile.write(str(eigenkets) + "\n")
    return None

# Decompose eigenvalList (with format [[k, eigenvals]]) into bands and make
# band plot.
def plotEigenvals(eigenvalList):
    if len(eigenvalList) == 0 or len(eigenvalList[0]) != 2:
        print("error: invalid input to plotEigenvals")
        return
    kpoints = []
    bands = []
    for i in range(len(eigenvalList[0][1])):
        bands.append([])
    for k, eigenvals in eigenvalList:
        kpoints.append(k)
        for i in range(len(eigenvals)):
            bands[i].append(eigenvals[i])
    for i in range(len(bands)):
        plt.plot(kpoints, bands[i])
    plt.show()

def main():
    # command line arguments
    calcType, kpointsFileName, outFileName = parseArgs()
    # read input file
    kpoints = getKpoints(kpointsFileName)
    # get appropriate Hamiltonian
    H = HamiltonianFn(calcType)

    seenZero = False # TODO - fix this hack - keeping only k_x for plot
    eigenvalList = []
    with open(outFileName, 'w') as outFile:
        # iterate over kpoints (parallelize later if necessary)
        for k in kpoints:
            # get Hamiltonian for this k-point
            Hk = H(k)
            # diagonalize Hamiltonian
            eigenvals, eigenkets = linalg.eigh(Hk)
            # handle output
            writeOutput(k, eigenvals, eigenkets, outFile)
            if k[1] == 0.0 and (k[0] != 0.0 or not seenZero): # TODO - fix this hack - keeping only k_x for plot
                if k[0] == 0.0:
                    seenZero = True
                a_hex = 4.138
                k[0] = 2.0 * math.pi * k[0] / a_hex
                eigenvalList.append([k[0], eigenvals])
    # plotting
    plotEigenvals(eigenvalList)

if __name__ == "__main__":
    main()
