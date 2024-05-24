all: sliding bending

sliding: Vars_sliding.f90 Mods_sliding.f90 Zconstrict_sliding.f90
	gfortran -fopenmp -fno-range-check $^ -o $@

bending: Vars_bending.f90 Mods_bending.f90 Zconstrict_bending.f90
	gfortran -fopenmp -fno-range-check $^ -o $@

