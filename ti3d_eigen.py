import sys
import json
import numpy as np
from scipy import linalg

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
    with open(calcType, 'r') as propsFile:
        props = json.load(propsFile)
        return fns[calcType](props)

# Generate 8-band Hamiltonian function with properties given by p
def Hamiltonian_8band(p):
    def H(k):
        return np.array([[1, 0], [0, -1]])
    return H

# Generate 4-band Hamiltonian function with properties given by p
def Hamiltonian_4band(p):
    def H(k):
        return np.array([[1, 0], [0, -1]])
    return H

# Generate mnk12 Hamiltonian function with properties given by p
def Hamiltonian_mnk12(p):
    def H(k):
        return np.array([[1, 0], [0, -1]])
    return H

# Write the output for one kpoint
def writeOutput(k, eigenvals, eigenkets, outFile):
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
