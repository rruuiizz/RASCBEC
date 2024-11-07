# RASCBEC
RAman Spectroscopy Calculation via Born Effective Charge

Code/rotate.py
POSCAR Rotation and Transformation Script

This script reads a VASP POSCAR file, applies specific rotational transformations to the atomic positions and lattice vectors, and outputs modified POSCAR files for three directions of applied electric fields: along the x, y, and z axes.

Input

POSCAR file: A VASP-formatted file named POSCAR which contains:
Lattice vector information.
Atom types and their quantities.
Atomic positions.
Output

Three VASP POSCAR-format files:

ex.POSCAR.vasp: Adjusted for electric field along the x-axis.
ey.POSCAR.vasp: Adjusted for electric field along the y-axis.
ez.POSCAR.vasp: Adjusted for electric field along the z-axis.
Each file includes the rotated lattice vectors and atomic positions.

