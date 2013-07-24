import sys, json, math
import numpy as np
from scipy import linalg
import matplotlib.pyplot as plt
import ti3d_eigen
from ti3d_eigen import Gamma1, Gamma2, Gamma3, Gamma4, Gamma5

# Parse command-line arguments numLayers, kpointsFileName, outFileName
def parseArgs():
    if len(sys.argv) != 4:
        print("Usage: python mnk12.py numLayers kpointsFileName outFileName")
        sys.exit(2)
    numLayers = int(sys.argv[1])
    if numLayers <= 0:
        print("Usage: python mnk12.py numLayers kpointsFileName outFileName")
        print("numLayers must be at least 1")
        sys.exit(2)

    return numLayers, sys.argv[2], sys.argv[3]

# Return a Hamiltonian function (k -> H_k) of the type specified
def HamiltonianFn(calcType, numLayers):
    fns = {"mnk12": Hamiltonian_mnk12}
    if calcType not in fns:
        print("error: calcType should be mnk12")
        sys.exit(2)
    with open(calcType + ".json", 'r') as propsFile:
        props = json.load(propsFile)
        return fns[calcType](props, numLayers)

# Generate mnk12 Hamiltonian function with properties given by p and the
# given number of layers.
def Hamiltonian_mnk12(p, numLayers):
    d = lambda k: p["M"] - 2.0*p["B"] + 2.0*p["B"] * (math.cos(k[0])
            + math.cos(k[1]) - 2.0)
    def H(k):
        # Convert k from 2pi/(lattice vector) units to 1/(lattice vector).
        # Lattice vector factor eliminated in sin(kx * a), etc.
        k = 2.0 * math.pi * k

        # 4x4 block along diagonal of 4Nx4N Hamiltonian (N=numer of layers)
        diagonal = np.diag([p["C"]]*4) + d(k)*Gamma5 + p["A"]*(Gamma2*math.sin(k[0]) + Gamma1*math.sin(k[1]))
        # 4x4 block appearing to the right of diagonal part (cross term)
        cross = p["B"]*Gamma5 - (1j*p["A"]/2.0)*Gamma4
        cross_conj = cross.T.conj()
        # 4Nx4N Hamiltonian
        Hk = np.zeros([4*numLayers, 4*numLayers], dtype=np.complex128)
        for i in range(0, 4*numLayers, 4): # i = 0, 4, 8, ..., 4*numLayers - 4
            Hk[i:i+4, i:i+4] = diagonal
            if i < 4*numLayers - 4: # exclude bottom-left cross
                Hk[i:i+4, i+4:i+8] = cross
            if i > 0:   # exclude top-right cross-conj
                Hk[i:i+4, i-4:i] = cross_conj
        return Hk

    return H

# Write the per-layer contribution of the given eigenket
def layerContribution(ket):
    layers = []
    for i in range(0, len(ket), 4):
        mag = 0.0
        for j in range(4):
            a = abs(ket[i+j])
            mag += a*a
        layers.append(mag)
    return layers

# Write the output for one kpoint
def writeOutput(k, eigenvals, eigenkets, outFile):
    #TODO - is there a better format?
    outFile.write("k:\n" + str(k) + "\n")
    outFile.write("eigenvals:\n" + str(eigenvals) + "\n")
    outFile.write("eigenkets:\n" + str(eigenkets) + "\n")
    outFile.write("layer contributions:\n")
    for ket in eigenkets:
        layers = layerContribution(ket)
        outFile.write(str(layers) + "\n")

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
    # get command line arguments
    numLayers, kpointsFileName, outFileName = parseArgs()
    # read input file
    kpoints = ti3d_eigen.getKpoints(kpointsFileName)
    # get appropriate Hamiltonian
    H = HamiltonianFn("mnk12", numLayers)

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
                eigenvalList.append([k[0], eigenvals])
    # plotting
    plotEigenvals(eigenvalList)

if __name__ == "__main__":
    main()
