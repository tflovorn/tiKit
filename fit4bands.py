import math
from numpy.linalg import norm
from scipy.optimize import leastsq
from scipy.linalg import eigh
from eigenval2foo import EIGENVAL
import ti3d_eigen

# Get the error corresponding to the params `p` at the specified k-point
# with the given 4-elem list of energy eigenvalues.
def H_err(p, k, energy):
    return norm(H_eigen(p, k) - energy)

# Get the expected 4-band eigenvalues at k-point `k` with params `p`.
def H_eigen(p, k):
    # convert list p to parameter map
    pmap = get_pmap(p)
    # get eigenvals for H_p(k)
    H_p = ti3d_eigen.Hamiltonian_4band(p)
    eigenvals, eigenkets = eigh(H_p(k))
    # convert eigenvals to numpy array
    return np.array(eigenvals)

# Convert `p`, a list of parameters, to a map representing these params.
def get_pmap(p):
    pmap = {"C0": p[0], "C1": p[1], "C2": p[2], "M0": p[3], "M1": p[4],
            "M2": p[5], "A0": p[6], "A2": p[7], "B0": p[8], "B2": p[9],
            "R1": p[10], "R2": p[11]}
    return pmap

# Get estimated list of parameters
def get_p_est():
    pass

# Get energies corresponding to the bands in the 4-band model from
# eigenvalue data. Assume these bands are the 2 above and 2 below the
# Fermi energy (verified that this assumption works with band data).
def getEnergyList(points):
    # get Fermi energy from OUTCAR

    # find 2 bands above and 2 bands below E_F

    # assemble list of eigenvalues for 4-band model
    pass

def main():
    # get eigenvalue data
    e = EIGENVAL("EIGENVAL")
    kList = e.kpoints
    energyList = getEnergyList(e.points)
    # get estimated parameters
    p_est = get_p_est()

    # perform least-squares fit for parameters
    p, ier = leastsq(H_err, p_est, args=(kList, energyList))
    if ier not in [1, 2, 3, 4]:
        # TODO handle error
        print("error in least-squares fit")

    # write output
    print(p)

if __name__ == "__main__":
    main()
