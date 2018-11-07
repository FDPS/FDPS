#!/bin/bash -x
# Check # of arguments
if [ $# -ne 1 ]; then
  echo "$# arguments are specified." 1>&2
  echo "This command requires single argument that represents the first ID of star particles." 1>&2
  exit 1
fi
# Set parameters
ID_STAR_1ST=$1
STAR_FILE_PREFIX="pos_star"
STAR_FILE_SUFFIX="txt"
SPH_FILE_PREFIX="pos_sph"
SPH_FILE_SUFFIX="txt"
regexp="([0-9]+)"
# Get file list
NBODY_FILES="nbody*.txt"
SPH_FILES="sph*.txt"
# Extract star particles from nbody*.txt
for fin in ${NBODY_FILES}; do
   if [[ ${fin} =~ ${regexp} ]];
   then
      fout="${STAR_FILE_PREFIX}${BASH_REMATCH[0]}.${STAR_FILE_SUFFIX}"
      awk -v id_1st=${ID_STAR_1ST} 'NR>2 {if ($1 >= id_1st) print $3, $4, $5}' < ${fin} > ${fout}
   fi
done
# Reduce SPH particle data
for fin in ${SPH_FILES}; do
   if [[ ${fin} =~ ${regexp} ]];
   then
      fout="${SPH_FILE_PREFIX}${BASH_REMATCH[0]}.${SPH_FILE_SUFFIX}"
      awk 'NR>2 {print $3, $4, $5}' < ${fin} > ${fout}
   fi
done
