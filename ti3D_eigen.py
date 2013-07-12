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
    #TODO open file or error
    kpoints = []
    #TODO read file (append to kpoints) or error
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
