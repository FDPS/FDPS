module physical_constants
   implicit none

   ! Physical constants
   double precision, parameter :: Ggrav = 6.67408d-8  ! Gravitational constants [cm^{3}/g/s^{2}]
   double precision, parameter :: kBoltz = 1.3806503d-16 ! Boltzmann constant [cm^{2} g/s^{2}/K]

   ! Mass units
   double precision, parameter :: Msolar = 1.989d33 ! Solar mass [g]
   double precision, parameter :: dalton = 1.660538921d-24 ! Unified atomic mass unit [g]
   double precision, parameter :: Mhydrogen = 1.007825d0 * dalton; ! Mass of hydrogen atom [g]

   ! Length units
   double precision, parameter :: km = 1.0d5 ! kilo-meters [cm]
   double precision, parameter :: AU = 1.49597870d13 ! Astronomical unit [cm]
   double precision, parameter :: pc = 3.0857d18 ! parsec [cm]
   double precision, parameter :: kpc = 1.0d3 * pc ! kilo-parsecs [cm]
   double precision, parameter :: Mpc = 1.0d6 * pc ! mega-parsecs [cm]
   double precision, parameter :: Gpc = 1.0d9 * pc ! giga-parsecs [cm]

   ! Time units
   double precision, parameter :: yr = 3.15576d7   ! year [s]
   double precision, parameter :: kyr = 1.0d3 * yr ! kilo-years [s]
   double precision, parameter :: Myr = 1.0d6 * yr ! mega-years [s]
   double precision, parameter :: Gyr = 1.0d9 * yr ! giga-years [s]

end module physical_constants
