if test $# -ne 2
then
    echo "$0 <nproc> <ifile>"
    exit
fi

nproc=$1
ifile=$2

export OMP_NUM_THREADS=1
mpirun-openmpi-gcc49 -np $nproc ./run $ifile >& log
theta=`grep tree_theta log | awk '{print $2}'`

for tag2 in pos app apm
do
    tag=theta"$theta"."$tag2"
    if test -e fdps_n32."$tag".0000
    then
        cat fdps_n32."$tag".* \
            | sort -n \
            > fdps_n32."$tag"
        rm -f fdps_n32."$tag".*
    fi
done

