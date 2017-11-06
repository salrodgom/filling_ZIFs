install:
	gfortran cif2lammps.f90 -o cif2lammps
	gfortran lammpstrj2pdb.f90  -o lammpstrj2pdb
	gfortran pdb2cif.f90  -o pdb2cif
	cc adsorption_fast_atom_saturation_01.c -lm -o adsorption_fast_atom_saturation_01
	cc zif_dist_angle_02.c -lm -o zif_dist_angle_02
clean:
	rm cif2lammps lammpstrj2pdb pdb2cif adsorption_fast_atom_saturation_01 zif_dist_angle_02
