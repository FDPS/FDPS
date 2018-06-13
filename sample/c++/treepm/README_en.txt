* How to Compile

Set the variables "FDPS_ROOT", "PM_ROOT", "FFTW_ROOT" in Makefile
depending on your computational envirnoment, perform "make" in the command-line.
If you want to use Phantom-GRAPE library, set the variable "PHANTOM_ROOT"
appropriately and make sure the variable "use_phantom_grape_x86" has the value "yes".


* How to use

Create a parameter file in the format described later, and perform
the executable "treepm" with specifying the absolute PATH of the 
parameter file as the 1st argument.

e.g.  mpiexec -n 4 treepm {the absolute PATH of the parameter file}

* Output files

Assuming the name of a run (model name) is "foo", this sample code
outputs the following files:
(1) "foo.diag", which records a simple runtime log.
(2) "foo_{output number}_{MPI rank number}" in the directory "foo_{output number}",
    which store snapshot of particle data for each MPI rank.
(3) "map_2d_{output number}", which stores the column number density of particles
    obtained by integrating the spacial number density over the z coordinate.


* Supported Formats of parameter file

This TreePM sample code supports the following three types of formats.


a) Format for reading an initial condition of Santa Barbara Cluster Comparison Test
/////////////////
0
Softening length for gravity
Opening angle criterion for the Tree method
Redshift at the end of a simulation
Name of snapshot files
Model name (name of a run)
Number of outputs of snapshots
Time of 1st snapshot Time of 2nd snapshot Time of 3rd snapshot ....
/////////////////
where the softening length is in the unit normalized by the side of
the computational box, the output times of snapshots should be given
as redshift values.

b) Format for reading a snapshot file of Santa Barbara Cluster Comparison Test
/////////////////
1
Softening length for gravity
Opening angle criterion for the Tree method
Redshift at the end of a simulation
Name of snapshot files
Redshift at the start of a simulation
\Omega_{M}
\Omega_{\Lambda}
\Omega_{b}
\Omega_{\nu}
Hubble constant in km/s/Mpc
Physical boxsize in Mpc/h where h is the Hubble parameter
Model name (name of a run)
Number of outputs of snapshots
Time of 1st snapshot Time of 2nd snapshot Time of 3rd snapshot ....
/////////////////
where the softening length is in the unit normalized by the side of
the computational box, the output times of snapshots should be given
as redshift values.

c) Format for an initial condition of randomly-distributed particles
/////////////////
2
Softening length for gravity
Opening angle criterion for the Tree method
Redshift at the end of a simulation
Redshift at the start of a simulation
\Omega_{M}
\Omega_{\Lambda}
\Omega_{b}
\Omega_{\nu}
Number of particles
Model name (name of a run)
Number of outputs of snapshots
Time of 1st snapshot Time of 2nd snapshot Time of 3rd snapshot ....
/////////////////
where the softening length is in the unit normalized by the side of
the computational box, the output times of snapshots should be given
as redshift values.
