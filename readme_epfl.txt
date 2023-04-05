need version: hdf5/1.10.4-intel19.1

some numerical inputs (#cpus,#grid points) defined in globals.f90

most input parameters defined in iTearing

the specifics/setup of the physics problem is defined in prob.f90

example of execution command: mpirun -n 4 ./hmhd2d Tearing

example of reading output from matlab: 
psi=h5read('data0004.hdf','/psi');
x=h5read('data0004,hdf','/x');
z=h5read('data0004,hdf','/z');
pcolor(z,x,psi');
