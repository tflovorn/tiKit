import sys
import json
import numpy as np
from scipy import linalg

# 3D model notation for Gamma matrices
Gamma1 = np.array([[0, 0, 0, 1],
                   [0, 0, 1, 0],
                   [0, 1, 0, 0],
                   [1, 0, 0, 0]])
Gamma2 = np.array([[0, 0, 0, -1j],
                   [0, 0, 1j, 0],
                   [0, -1j, 0, 0],
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
    # iterate over (kStart, kStop) pairs
    for pair in kBounds:
        start, stop, step = pair[0], pair[1], []
        # possibly different step for each coordinate
        for i in range(len(start)):
            step.append((stop[i] - start[i])/(kNum-1))
        # make kNum points between start and stop, including boundaries
        for ptIndex in range(kNum):
            point = []
            for i in range(len(start)):
                point.append(start[i] + ptIndex*step[i])
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
    #TODO
    print("Warning: 8band Hamiltonian is unimplemented")
    def H(k):
        return np.array([[1, 0], [0, -1]])
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
    A = lambda k: p["A0"] - 0.5*p["A0"]*kp2(k) #TODO check guesses A2 = (-1/2) A0 and B2 = (-1/2) B0
    B = lambda k: p["B0"] - 0.5*p["B0"]*kz2(k)
    # parenthetic expression in H3, Eq. 17, Liu 2010
    q = lambda kx, ky: kx**3 - 3.0*kx*(ky**2)

    def H(k):
        H0 = epsilon(k) + M(k)*Gamma5 + B(k)*Gamma4*k[2] + A(k)*(Gamma1*k[1] - Gamma2*k[0])
        H3 = p["R1"]*Gamma3*q(k[0], k[1]) - p["R2"]*Gamma4*q(k[1], k[0])
        return H0 + H3

    return H

# Generate mnk12 Hamiltonian function with properties given by p
def Hamiltonian_mnk12(p):
    #TODO
    print("Warning: mnk12 Hamiltonian is unimplemented")
    def H(k):
        return np.array([[1, 0], [0, -1]])
    return H

# Write the output for one kpoint
def writeOutput(k, eigenvals, eigenkets, outFile):
    #TODO
    print(k)
    print(eigenvals)
    print(eigenkets)
    return None

def main():
    # command line arguments
    calcType, kpointsFileName, outFileName = parseArgs()
    # read input file
    kpoints = getKpoints(kpointsFileName)
    # get appropriate Hamiltonian
    H = HamiltonianFn(calcType)

    with open(outFileName, 'w') as outFile:
        # iterate over kpoints (parallelize later)
        for k in kpoints:
            Hk = H(k)
            eigenvals, eigenkets = linalg.eig(Hk)
            writeOutput(k, eigenvals, eigenkets, outFile)

if __name__ == "__main__":
    main()
