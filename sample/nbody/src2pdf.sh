#!/bin/sh

awk '{printf("%-4d %s\n", NR, $0)}' user-defined.hpp > user-defined.txt
a2ps --no-header --center-title=user-defined.hpp --pretty-print=cxx user-defined.txt -o user-defined.ps
ps2pdf user-defined.ps
rm user-defined.ps user-defined.txt
awk '{printf("%-4d %s\n", NR, $0)}' nbody.cpp > nbody.txt
a2ps --no-header --center-title=nbody.cpp --pretty-print=cxx nbody.txt -o nbody.ps
ps2pdf nbody.ps
rm nbody.ps nbody.txt
