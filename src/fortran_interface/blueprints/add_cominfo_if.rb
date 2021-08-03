#
# add_cominfo_if_rb
#
# create new  (Version 7)  FDPS_c_if_blueprint.h
# from old one (version 6)
#
in_comm = false
header_printed = false
while s=gets
  in_comm = true  if s =~ /MPI comm/
  in_comm = false  if s =~ /Utility/
  unless in_comm
    print s.chomp, "\n"
  else
    a=s.split
    if a[1] =~/^fdps_/
      while s.index(";")==nil
        s+=gets
      end
      snew=s.gsub(/(\S+)\s+(fdps_)(\S+\()/){
        $1+" fdps_ci_"+$3+"fdps_comm_info * ci, "}
      snew.gsub!(/,\s*\)/,")")
      print snew
      snew=s.gsub(/(\S+)\s+(fdps_)(\S+\()/){
        $1+" fdps_f_ci_"+$3+"int ci, "}
      snew.gsub!(/,\s*\)/,")")
      print snew
    end
    print s
    if s =~ /----------/ && !header_printed
      header_printed = true
      print <<-EOF

#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
#include <mpi.h>
#endif
#include "FDPS_comm_info.h"

#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
#define MPI_Comm int
#define MPI_Fint int
#endif    
fdps_comm_info * fdps_ci_initialize(MPI_Comm* comm);
void fdps_ci_set_communicator(fdps_comm_info * ci, MPI_Comm * comm);

int fdps_f_ci_initialize(MPI_Fint comm);
void fdps_f_ci_set_communicator(int ci, MPI_Fint comm);
void fdps_f_ci_delete(int ci);




void fdps_ci_delete(fdps_comm_info * ci);
fdps_comm_info * fdps_ci_create(fdps_comm_info * ci, int n, int rank[]);
fdps_comm_info * fdps_ci_split(fdps_comm_info * ci,int color, int key);

int  fdps_f_ci_create(int ci, int n, int rank[]);
int fdps_f_ci_split(int ci,int color, int key);

EOF
    end  
  end
end


    
