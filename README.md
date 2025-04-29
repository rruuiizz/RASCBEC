# RASCBEC

RAman Spectroscopy Calculation via Born Effective Charge

# Overview of RASCBEC Method

This method calculates Raman activities based on Born Effective Charge (BEC) data, following phonon calculations. Two separate scripts are provided depending on how the phonon modes were obtained:

- RASCBEC_vasp.py — for phonon data generated directly by VASP
- RASCBEC_phonopy.py — for phonon data generated using phonopy

# Workflow Steps

- Rotation of Structures
Use the script rotate.py to generate rotated POSCAR files along ±x, ±y, and ±z directions.
This prepares the necessary configurations for BEC calculations.

- Born Effective Charge (BEC) Calculations
Perform 8 separate VASP calculations using the rotated structures.
Extract the BEC tensors from each OUTCAR file.

- Raman Activity Computation
Run RASCBEC_vasp.py if phonons were computed via VASP.
Run RASCBEC_phonopy.py if phonons were computed via phonopy.
These scripts will combine phonon frequencies and eigenvectors with the BEC derivatives to compute Raman activities.

# Detailed Description of Each Script

- rotate.py (POSCAR Rotation and Transformation Script)

This script reads a VASP POSCAR file, applies specific rotational transformations to the atomic positions and lattice vectors, and outputs modified POSCAR files for three directions of applied electric fields: along the x, y, and z axes.

- Input

POSCAR file: A VASP-formatted file named POSCAR which contains:
Lattice vector information.
Atom types and their quantities.
Atomic positions.

- Output

Three VASP POSCAR-format files:
ex.POSCAR.vasp: Adjusted for electric field along the x-axis.
ey.POSCAR.vasp: Adjusted for electric field along the y-axis.
ez.POSCAR.vasp: Adjusted for electric field along the z-axis.
Each file includes the rotated lattice vectors and atomic positions.

- RASCBEC_vasp.py (Calculation of Raman Activities from VASP Output Files)
 
This script reads input files generated from VASP simulations 
to calculate Raman activities for each phonon mode using RASCBEC method.

It requires:

- A POSCAR file (structure information).
- Eight OUTCAR files (OUTCAR1, OUTCARm1, OUTCARx, OUTCARmx, OUTCARy, OUTCARmy, OUTCARz, OUTCARmz) containing Born Effective Charge data.
- Two additional phonon property files:
    - freqs_vasp.dat: Phonon frequencies stored as a 3N×1 array, where N is the number of atoms.
    - eigvecs_vasp.dat: Phonon eigenvectors stored as a 3N×3N array.

Dependencies:
- numpy

Usage:

Prepare the following files in the same directory:
- POSCAR
- OUTCAR1, OUTCARm1, OUTCARx, OUTCARmx, OUTCARy, OUTCARmy, OUTCARz, OUTCARmz
- freqs_vasp.dat
- eigvecs_vasp.dat
Run the script using Python 3 and NumPy installed:
python RASCBEC_vasp.py

Output:

- - raman_vasp.dat:
A file containing the calculated Raman activities for each phonon mode.

- RASCBEC_vasp.py
Similar to RASCBEC_phonopy but requiring freqs_phonopy.dat and eigvecs_phonopy.dat


# Example Case: Rutile GeO2

In this example, we demonstrate the full Raman activity calculation workflow using rutile GeO2, which contains 6 atoms per unit cell.

All necessary input and output files are provided in the example folder:

Included Files

Input Files:
POSCAR — structure of rutile GeO2
OUTCAR1, OUTCARm1, OUTCARx, OUTCARmx, OUTCARy, OUTCARmy, OUTCARz, OUTCARmz — BEC calculations for 8 rotated structures
freqs_vasp.dat — phonon frequencies in a 18×1 array format (N = 6) got from VASP
eigvecs_vasp.dat — phonon eigenvectors in a 18×18 array got from VASP
freqs_phonopy.dat — phonon frequencies in a 18×1 array format (N = 6) got from phonopy
eigvecs_phonopy.dat — phonon eigenvectors in a 18×18 array got from phonopy

Output File:
raman_vasp.dat — computed Raman activities for each phonon mode for phonon data generated directly by VASP
raman_phonopy.dat — computed Raman activities for each phonon mode for phonon data generated using phonopy

This example can be used as a reference for setting up and validating your own Raman activity calculations.


# Citation: 
Zhang, Rui, et al. "RASCBEC: Raman spectroscopy calculation via born effective charge." Computer Physics Communications 307 (2025): 109425

