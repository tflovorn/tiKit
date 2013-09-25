Resources for topological insulator calculations, including ab initio inputs and model code.

All content released under the MIT license except for vasp\_scripts/eigenval2foo.py, which predates this project and is released under the GPLv2.

Ab initio inputs:

* vasp\_scripts/bi2se3poscar: generates POSCAR for Bi2Se3 hexagonal or rhombohedral unit cell (number of layers is a multiple of 15).

* vasp\_scripts/bi2se3poscar\_frac: generates POSCAR for Bi2Se3 hexagonal cell (number of layers is a multiple of 5).

* vasp\_scripts/bi2se3poscar\_surf\_relaxed: generates POSCAR for Bi2Se3 hexagonal slab using published data for relaxed distance between quintuple-layers.

* vasp\_inputs/Bi2Se3/bulk: INCAR, KPOINTS, and POSCAR files for bulk rhombohedral and hexagonal cell band structure.

    * The POTCAR file for Bi and Se must be provided. We use the PAW\_PBE POTCAR files provided with VASP.

    * Subdirectories are organized as \*-chg and \*-bands. The chg inputs generate the charge density files CHG, CHGCAR which are then copied into the bands directory.

Model code:

* ti3D\_eigen/ti3d\_eigen.py & ti3D\_eigen/mnk12.py: implement bulk and slab Hamiltonians for Bi2Se3. For more details see ti3D\_eigen/README.md.

* ti3D\_gf: Fortran implementation of slab Hamiltonian. Will be used to obtain Green's functions for performing transport calculations.
