Implements diagonalization of the Hamiltonians modeling topological insulators from the following papers:

[Liu et al PRB 82, 045122 (2010)](http://prb.aps.org/abstract/PRB/v82/i4/e045122)

[Mahfouzi, Nagaosa and Nikolic PRL 109, 166602 (2012)](http://prl.aps.org/abstract/PRL/v109/i16/e166602)

Dependencies
==============

Requires scipy for matrix diagonalization and matplotlib for plotting. In Debian and its derivatives these can be obtained with:

    sudo apt-get install python-scipy python-matplotlib python-tk

Usage
==========

The bulk diagonalization can be run with the command:

    python ti3d_eigen.py TYPE KPOINTS OUTFILE

Here TYPE = 4band, 8band, or mnk12; KPOINTS = path of the input file with k-points to diagonalize H over (VASP band calculation format - see KPOINTS\_symmetry for an example); OUTFILE = path of the file to write output to.
The 4band and 8band Hamiltonians are in the continuum limit and taken from Liu (2010). The mnk12 Hamiltonian is tight-binding and is taken from Mahfouzi (2012).

The diagonalization of the system with a finite number of layers can be run with the command:

    python mnk12.py NUM_LAYERS KPOINTS OUTFILE

Here NUM\_LAYERS is the number of layers in the finite-layer system and KPOINTS and OUTFILE are as before.
