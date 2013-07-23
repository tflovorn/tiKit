Implements diagonalization of the Hamiltonians modeling topological insulators from the following papers:

[Liu et al PRB 82, 045122 (2010)](http://prb.aps.org/abstract/PRB/v82/i4/e045122)

[Mahfouzi, Nagaosa and Nikolic PRL 109, 166602 (2012)](http://prl.aps.org/abstract/PRL/v109/i16/e166602)

Usage
==========

Requires scipy to be installed (package python-scipy in Debian derivatives).

Run with the command:

    python ti3d_eigen.py TYPE KPOINTS OUTFILE

Here TYPE = 4band, 8band, or mnk12; KPOINTS = path of the input file with k-points to diagonalize H over (VASP band calculation format - see KPOINTS\_symmetry for an example); and OUTFILE = path of the file to write output to.

The 4band and 8band Hamiltonians are in the continuum limit and taken from Liu (2010). The mnk12 Hamiltonian is tight-binding and is taken from Mahfouzi (2012).
