#
# add_cominfo_if_cpp.rb
#
# create new  (Version 7)  FDPS_c_if_blueprint.cpp
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
      ss=s
      while ss.index("}")!= 0
        ss=gets
        s+=ss
      end
      snew=s.gsub(/(\S+)\s+(fdps_)(\S+\()/){
        $1+" fdps_ci_"+$3+"fdps_comm_info * ci, "}
      snew.gsub!(/,\s*\)/,")")
      snew.gsub!(/PS::Comm::/,"((PS::CommInfo*)(ci))->");
      print snew
      snew=s.gsub(/(\S+)\s+(fdps_)(\S+\()/){
        $1+" fdps_f_ci_"+$3+"int ci, "}
      snew.gsub!(/,\s*\)/,")")
      snew.gsub!(/PS::Comm::/,"((PS::CommInfo*)(comm_table[ci]))->");
      print snew
    end
    print s
    if s =~ /----------/ && !header_printed
      header_printed = true
      print <<-EOF
#include "FDPS_comm_info.h"


void fdps_ci_delete(fdps_comm_info * ci)
{
    delete ((PS::CommInfo*)ci);
}

const int fdps_comm_table_size = 32;
static int fdps_comm_table_initialized = 0;
static fdps_comm_info * comm_table[fdps_comm_table_size];

void fdps_ci_initialized_comm_table_if_not()
{
   if (fdps_comm_table_initialized ==0){        
       for(int i=0;i<fdps_comm_table_size; i++){
           comm_table[i]=NULL;
       }
       fdps_comm_table_initialized = 1;
    }
}

int  fdps_ci_register_comm(fdps_comm_info * p)
{
   fdps_ci_initialized_comm_table_if_not();
   for(int i=0;i<fdps_comm_table_size; i++){
       if (comm_table[i]==NULL){
           comm_table[i]==p;
           return i;
       }
   }
   return -1;
}

#ifndef PARTICLE_SIMULATOR_MPI_PARALLEL
#define MPI_Comm int
#define MPI_Fint int

#endif
fdps_comm_info * fdps_ci_initialize(MPI_Comm* comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    PS::CommInfo CI = PS::CommInfo( *comm);
#else
    PS::CommInfo CI = PS::CommInfo();
#endif
    return (fdps_comm_info *) (&CI);
}
void fdps_ci_set_communicator(fdps_comm_info * ci,MPI_Comm * comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    ((PS::CommInfo*)ci)->setCommunicator(*comm);
#else
    ((PS::CommInfo*)ci)->setCommunicator();
#endif

}


int fdps_f_ci_initialize(MPI_Fint comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Comm ccomm = MPI_Comm_f2c(comm);
#else
    MPI_Comm ccomm = comm;
#endif
    return fdps_ci_register_comm(fdps_ci_initialize(&ccomm));
}
void fdps_f_ci_set_communicator(int ci, MPI_Fint comm)
{
#ifdef PARTICLE_SIMULATOR_MPI_PARALLEL
    MPI_Comm ccomm = MPI_Comm_f2c(comm);
#else
    MPI_Comm ccomm = comm;
#endif
    fdps_ci_set_communicator(comm_table[ci], &ccomm);
}

void fdps_f_ci_delete(int ci)
{
    fdps_ci_delete(comm_table[ci]);
    comm_table[ci] = NULL;
}
fdps_comm_info * fdps_f_ci_getcomtable(int i)
{
   return (fdps_comm_info * ) comm_table[i];
}



fdps_comm_info * fdps_ci_create(fdps_comm_info * ci, int n, int rank[])
{
    PS::CommInfo CI =((PS::CommInfo*)ci)->create(n, rank);
    return (fdps_comm_info *) (&CI);
}
fdps_comm_info * fdps_ci_split(fdps_comm_info * ci,int color, int key)
{
    PS::CommInfo CI =((PS::CommInfo*)ci)->split(color, key);
    return (fdps_comm_info *) (&CI);
}
int  fdps_f_ci_create(int ci, int n, int rank[])
{
    return fdps_ci_register_comm(fdps_ci_create(comm_table[ci], n, rank));
}
int fdps_f_ci_split(int ci,int color, int key)
{
    return fdps_ci_register_comm(fdps_ci_split(comm_table[ci], color, key));
}



EOF
    end  
  end
end


    
