#!/bin/bash

#-DNELE=$power
CFLAGS="-O3 -funroll-loops -floop-unroll-and-jam -march=cascadelake -mavx512f -mavx512dq -mavx512ifma"

for p in {30..6..-1}
do 
    power=$((2 ** $p))
    g++ $CFLAGS -DTYPE=double -DNSAMPLE=20 -DOPTIM=transform -o fmadd fmadd.cpp
    ./fmadd $power
done

for p in {30..6..-1}
do 
    power=$((2 ** $p))
    g++ $CFLAGS -DTYPE=double -DNSAMPLE=20 -DOPTIM=mipp -o mippfmadd mipp_fmadd.cpp
    ./mippfmadd $power
done

echo "DONE"