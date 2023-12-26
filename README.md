# PyMoire
A package for tight-binding calculation of twisted bilayer graphene based on mapped Wannier functions
The PyMoire package calculates the band structure and DOS for twisted bilayer graphene in any commensurate twist angle with the tight-binding model
using the Hopping parameters extracted from the Wannier functions based on the DFT calculation in the Quantum-Espresso package.
The Wannier hopping parameters are calculated in different stacking configurations and then matched on the moire structure in real space.

----------------------------------------------------
To obtain hopping parameters with respect to inter-atomic distance, we displace one layer relative to another layer in a specified steps in QE_dxdy_AA_displacing script file.

Please modify the QE_dxdy_AA_displacing file. Determine Quantum-Espresso/bin, wannier, and pseudo-potential directories. For AA stacking run the script in the Linux terminal:
sh ./QE_dxdy_AA_displacing.sh
and for AB stacking run:
sh ./QE_dxdy_AB_displacing.sh
If file needs to permission enter this command:
chmod +x QE_dxdy_AA_displacing.sh

In each step, scf, nscf and Wannier calculations are performed and the results are extracted and stored in files that include: Total Energy, VdW Energy, Fermi Energy, Band gap (the scf file should be modified) and hopping parameters of coupled orbitals individually. The hopping parameter files are named after paired orbital numbers which are determined in wannier input section. For example in this calculation 4-9.dat is pz of layer 1 to pz of layer 2.

The first and second columns of hopping parameter files are displacement of dx and dy form origin and the rest of the columns are hopping parameters in Wigner–Seitz cell order which we choose 0,0,0 cell.

This calculation takes several hours (depending on the number of parallel cores used). By default, the files required for tight-binding calculations are placed in the ./out folder (If you run the script, delete this folder) 

----------------------------------------------------
The PyMoire code is programmed for twisted bilayre graphene in python language which the main file is PyMoire.py and it uses functions in the PyMoireFunc.py file. The libraries used in this model are:
Numpy, Scipy, Matplotlib, libtetrabz (DOS calculation) and time (timing of each calculation section).

The model is structured as: 
1- Constructing bilayer graphene and twisting one layer respect to another and calculating moire lattice vectors, unit-cell, and Brillouin zone.
2- Calculate interlayer coupling distance for three nearest neighbors atoms.
3- Entering hopping parameters and removing outlier data with K-nearest-neighbor method.
4- Interpolate the hopping parameters with Radial Basis Function method and mapping them onto moire structure.
5- Indexing coupled orbitals.
6- Sorting and stacking orbital indexes and corresponding hopping and coupling vectors.
7- Constructing tight-binding Hamiltonian and solving it.
8- Exporting hr.dat, bandstructure and DOS files and plotting the results.

To calculate different twist angles of the default structure (twisted bilayer graphene), simply change the twist constant number in PyMoire file (described in the code). The code automatically finds the minimum number of atoms to assign to the moire supercell. Calculates lattice constant and Brillouin zone. Then calculates the interlayer distances of atoms and matches the hopping values based on the hopping maps calculated from DFT/Wannier and finally completes tight-binding calculations.

To calculate different structures, you need to run the DFT/Wannier calculation for the corresponding structure and then modify the structure section in the PyMoire file. You may need to change some functions as well. But the general calculation algorithm, from structure twisting to parameters mapping and tight-binding calculations, is suitable for all structures.

By default, the results of twisted bilayer graphene calculations for a twist angle of 3.48 and 1.05 degrees are placed in the ./out folder including plot of mapped moire sturcture, band structure and DOS data file and their plots, also, the hr.dat file that can be used in other packages like WannierTools and PythTB.
In the future, more details such as spin-orbit calculations and the Hubbard model will be added.
-----------------------------------------------------
Citation

Pymoire is an open-source Python package designed for the analysis of moiré patterns.
If you use the Pymoire package in your work, we kindly request you to cite the following reference:
https://doi.org/10.1016/j.physe.2023.115877
BibTex:
@article{servati2023real,
  title={Real-space tight-binding model for twisted bilayer graphene based on mapped Wannier functions},
  author={Servati, Mahyar and Rasuli, Reza and Tavana, Ali},
  journal={Physica E: Low-dimensional Systems and Nanostructures},
  pages={115877},
  year={2023},
  publisher={Elsevier}
}


