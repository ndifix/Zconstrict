# Zconstrict

To install a Fortran compiler in Ubuntu, run the following in a terminal:
sudo apt-get install gfortran

To compile ZCONSTRICT, run the following in a terminal (which takes ~ 1 min):

for the sliding model:

gfortran -fopenmp -fno-range-check Vars_sliding.f90 Mods_sliding.f90 Zconstrict_sliding.f90 -o runfile


for the bending model:

gfortran -fopenmp -fno-range-check Vars_bending.f90 Mods_bending.f90 Zconstrict_bending.f90 -o runfile



To run simulation:

./runfile >log.txt

Expected output: zring.psf, startzring.pdb, traj*.dcd

Note: the software has been tested on Linux operating systems including Ubuntu 16.04.4 LTS


The data can be visualized using VMD, which can be downloaded from: https://www.ks.uiuc.edu/Research/vmd/

in a terminal, run the following to visualize the output:
vmd *.psf and *.dcd

