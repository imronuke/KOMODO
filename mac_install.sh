#!/bin/bash
# start

echo "KOMODO compilation starts"

set FORT=gfotran

echo " "
echo "Compiling src/mod_data.f90"
gfortran -O4 -c src/mod_data.f90
echo "Compiling src/mod_io.f90"
gfortran -O4 -c src/mod_io.f90
echo "Compiling src/mod_xsec.f90"
gfortran -O4 -c src/mod_xsec.f90
echo "Compiling src/mod_nodal.f90"
gfortran -O4 -c src/mod_nodal.f90
echo "Compiling src/mod_cmfd.f90"
gfortran -O4 -c src/mod_cmfd.f90
echo "Compiling src/mod_th.f90"
gfortran -O4 -c src/mod_th.f90
echo "Compiling src/mod_trans.f90"
gfortran -O4 -c src/mod_trans.f90
echo "Compiling src/mod_control.f90"
gfortran -O4 -c src/mod_control.f90
echo "Compiling src/komodo.f90"
gfortran -O4 -c src/komodo.f90
echo "Combining all together"
gfortran *.o -o komodo -static-libgfortran

echo " "
echo "Copy komodo to /usr/local/bin"
sudo cp komodo /usr/local/bin

echo " "
echo "Deleted unnecessary files"
rm *.o *.mod

echo " "
echo "KOMODO successfully compiled"
