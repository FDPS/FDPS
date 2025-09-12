* How to Compile

Edit "FDPS_ROOT" and "FFTW_ROOT" in the Makefile as appropriate, then run make.
If you want to use pikg, set "use_pikg=yes" and set "PIKG_ROOT" as well.

* How to Use

Create a parameter file and specify its absolute PATH as the first argument.

 e.g.  mpiexec -n 4 pmmm.out {parameter file}

* Parameter File Format

This sample code supports the following three formats.
If the model name is "foo", a simple output log will be written to "foo.diag".
Snapshots will be output to files named "foo_{output number}-{process rank}" in the directory "foo_{output number}" (one file per process),
and the column density in the x-y plane will be output to a file named "map_2d_{output number}".

a) When reading the initial conditions of the Santa Barbara Cluster
/////////////////
0
Softening length between particles
Tree opening angle
Final z
Snapshot file name
Model name
Number of snapshot outputs
Snapshot time 0 Snapshot time 1 Snapshot time 2 ....
/////////////////
Here, the softening length is given in units where the box size is 1.
Snapshot times should be given as redshift values.

b) When reading a snapshot file in the same format as the Santa Barbara Cluster
/////////////////
1
Softening length between particles
Tree opening angle
Final z
Snapshot file name
Initial z
Omega matter
Omega lambda
Omega baryon
Omega radiation
Hubble constant (km/s/Mpc)
Box size (MPC/h)
Model name
Number of snapshot outputs
Snapshot time 0 Snapshot time 1 Snapshot time 2 ....
/////////////////
Here, the softening length is given in units where the box size is 1.
Snapshot times should be given as redshift values.

c) Simulation with randomly distributed particles
/////////////////
2
Softening length between particles
Tree opening angle
Final z
Initial z
Omega matter
Omega lambda
Omega baryon
Omega radiation
Number of particles
Model name
Number of snapshot outputs
Snapshot time 0 Snapshot time 1 Snapshot time 2 ....
/////////////////
Here, the softening length is given in units where the box size is 1.
Snapshot times should be given as redshift values.
