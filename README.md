# test-mipp
testing mipp

# How to run
## Preprocessor macros
Available at compile time:
- TYPE in double, float or int, the data type used
- NSAMPLE integer, running the kernel function to be bench NSAMPLE times
- OPTIM in transform, mipp, the kernel function to be used
## Input args
Available at run time:
- NELE integer, vector length
## Example
```bash
#!/bin/bash

guix shell gcc-toolchain@13.2.0
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
```
