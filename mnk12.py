import sys, json, math
import numpy as np
from scipy import linalg
import ti3d_eigen
from ti3d_eigen import Gamma1, Gamma2, Gamma3, Gamma4, Gamma5

# Parse command-line arguments numLayers, kpointsFileName, outFileName
def parseArgs():
    if len(sys.argv) != 4:
        print("Usage: python mnk12.py numLayers kpointsFileName outFileName")
        sys.exit(2)
    return int(sys.argv[1]), sys.argv[2], sys.argv[3]

# Return a Hamiltonian function (k -> H_k) of the type specified
def HamiltonianFn(calcType):
    fns = {"mnk12": Hamiltonian_mnk12}
    if calcType not in fns:
        print("error: calcType should be mnk12")
        sys.exit(2)
    with open(calcType + ".json", 'r') as propsFile:
        props = json.load(propsFile)
        return fns[calcType](props)

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

def main():
    # command line arguments
    numLayers, kpointsFileName, outFileName = parseArgs()
    # read input file
    kpoints = ti3d_eigen.getKpoints(kpointsFileName)
    # get appropriate Hamiltonian
    H = HamiltonianFn("mnk12")

    with open(outFileName, 'w') as outFile:
        # iterate over kpoints (parallelize later if necessary)
        for k in kpoints:
            Hk = H(k)
            eigenvals, eigenkets = linalg.eigh(Hk)
            writeOutput(k, eigenvals, eigenkets, outFile)

if __name__ == "__main__":
    main()
