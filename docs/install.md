---
title: Quick Install Guide
theme: _config.yml
filename: install
---

# Building from source codes
To build ADPRES from the source codes you just need a Fortran compiler, that's it. We design ADPRES to be very portable, so we expect it can be installed using any Fortran compiler on any machine.

## Building in Ubuntu or other GNU-Linux based OS
In Ubuntu, other GNU-Linux based OS or CYGWIN you can use either gfotran or Intel fortran to build the source codes. You can install gfortran in Ubuntu OS by using this command in your computer's terminal

```
sudo apt install gfortran
```

Then you can download the [ADPRES zip files](https://github.com/imronuke/ADPRES/archive/master.zip) then extract that zip file, or you can clone them (note: you need to download git first if you don't have it)

```
git clone https://github.com/imronuke/ADPRES.git
```

In a machine where the gfortran is already installed, go to the ADPRES folder which you had been extracted or cloned

```
cd ADPRES-master
```

and build the source codes by using command:

```
sudo ./install.sh
```

If this way didn't work, try to install bash into your computer. Alternatively, you can build ADPRES manually

```
gfortran -O4 -c src/mod_data.f90
gfortran -O4 -c src/mod_io.f90
gfortran -O4 -c src/mod_xsec.f90
gfortran -O4 -c src/mod_nodal.f90
gfortran -O4 -c src/mod_cmfd.f90
gfortran -O4 -c src/mod_th.f90
gfortran -O4 -c src/mod_trans.f90
gfortran -O4 -c src/mod_control.f90
gfortran -O4 -c src/ADPRES.f90
gfortran *.o -o adpres
sudo cp adpres /usr/bin
```

(NOTE: you need to have admin privilege to run these commands)

These command will create an executable file named `adpres` and copied the executable file to `/usr/bin`. Now, you can run a test using several examples of inputs file in folder [smpl](https://github.com/imronuke/ADPRES/tree/master/smpl) to see if you had built ADPRES properly. You can run ADPRES using command

```
adpres [INPUT_FILE_PATH_NAME]
```

for example, you can run [`IAEA3Ds`](https://github.com/imronuke/ADPRES/blob/master/smpl/static/IAEA3Ds) input by

```
adpres smpl/static/IAEA3Ds
```

If you see `ADPRES EXIT NORMALLY` at the end of terminal output, then congratulations! you have successfully installed ADPRES. By the way, it should take about 0.2 seconds for ADPRES to solve IAEA3Ds problem if you build using gfortran in a typical today's computer. The way to build using Intel fortran is similar, just change `gfortran` with `ifort` in the commands above.

## Building in Mac OS
The way to install ADPRES on Mac OS is similar to installing ADPRES on Ubuntu. First you need to install gfortran to your computer. You can find the information on how to install gfortran for Mac OS from [here](https://gcc.gnu.org/wiki/GFortranBinariesMacOS). And instead of executing this command as you were supposed to do in Ubuntu

```
sudo ./install.sh
```

you shall execute this command in Mac OS

```
sudo ./mac_install.sh
```

The rest of the steps are the same.


## Building in Windows
There are at least two ways installing ADPRES on Windows.

### Using Windows Subsystem for Linux (WSL) in Windows 10
Using this way, you can go to Microsoft Store and install Ubuntu from there for free. Open the Ubuntu app and follow the same steps on installing ADPRES on Ubuntu OS.

### Using g95 Fortran Compiler
You can also intsall ADPRES directly into Windows by using free fortran compiler g95 which can be obtained from [here](https://www.fortran.com/wp-content/uploads/2013/05/g95-Mingw_201210.exe). Then install g95 to your computer.

Now you can download the [ADPRES zip files](https://github.com/imronuke/ADPRES/archive/master.zip) from Github then extract that zip file.

After you installed g95 and extracted ADPRES zip file, open the command prompt. And from the command prompt, using `cd` command, go to the  ADPRES folder which you had been extracted

```
cd ADPRES-master
```

and build the source codes using g95:

```
g95 -O4 -c src\mod_data.f90
g95 -O4 -c src\mod_io.f90
g95 -O4 -c src\mod_xsec.f90
g95 -O4 -c src\mod_nodal.f90
g95 -O4 -c src\mod_cmfd.f90
g95 -O4 -c src\mod_th.f90
g95 -O4 -c src\mod_trans.f90
g95 -O4 -c src\mod_control.f90
g95 -O4 -c src\ADPRES.f90
g95 *.o -o adpres
```

These command will create an executable file named `adpres`. If you use Intel fortran compiler, just change `g95` with `ifort` in the commands above. Now, you can run a test using several examples of inputs file in folder [smpl](https://github.com/imronuke/ADPRES/tree/master/smpl) to see if you had built ADPRES properly. You can run ADPRES using command

```
adpres [INPUT_FILE_PATH_NAME]
```

for example, you can run [`IAEA3Ds`](https://github.com/imronuke/ADPRES/blob/master/smpl/static/IAEA3Ds) input by

```
adpres smpl\static\IAEA3Ds
```

If you see `ADPRES EXIT NORMALLY` at the end of terminal output, then congratulations! you have successfully installed ADPRES on Windows.
