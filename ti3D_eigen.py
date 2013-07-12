import sys
import numpy as np
from scipy import linalg

# Parse command-line arguments calcType, kpointsFileName, outFileName
def parseArgs():
    if len(sys.argv) != 4:
        print("Usage: ti3D_eigen [8band|4band|mnk12] kpointsFileName outFileName")
        sys.exit(2)
    return sys.argv[1], sys.argv[2], sys.argv[3]

# Parse kpoints file
def getKpoints(kpointsFileName):
    #TODO read file (append to kpoints) or error
    #kpointsFile = open(kpointsFileName, 'r')
    kNum = 25
    kBounds = [[[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]], [[0.5, 0.5, 0.5], [0.0, 0.5, 0.5]]]
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
    # valid calcTypes: "8band", "4band", "mnk12"
    # TODO implement (8band first)
    def H(k):
        return 0.0
    return H

# Return the eigenvalues and corresponding eigenvectors of M
def EigenDecompose(M):
    #TODO
    return [], []

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
            eigenvals, eigenkets = EigenDecompose(Hk)
            writeOutput(k, eigenvals, eigenkets, outFile)

if __name__ == "__main__":
    main()
