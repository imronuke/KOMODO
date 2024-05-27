---
title: Quick Install Guide
theme: _config.yml
filename: install
---

# Building from source codes
To build KOMODO from the source codes you just need a Fortran compiler, that's it. We design KOMODO to be very portable, so we expect it can be installed using any Fortran compiler on any machine.

## Building in Ubuntu or other GNU-Linux based OS
In Ubuntu, other GNU-Linux based OS or CYGWIN you can use either gfotran or Intel fortran to build the source codes. You can install gfortran in Ubuntu OS by using this command in your computer's terminal

```
sudo apt install gfortran
```

Then you can download the [KOMODO zip files](https://github.com/imronuke/KOMODO/archive/master.zip) then extract that zip file, or you can clone them (note: you need to download git first if you don't have it)

```
git clone https://github.com/imronuke/KOMODO.git
```

In a machine where the gfortran is already installed, go to the KOMODO folder which you had been extracted or cloned

```
cd KOMODO-master
```

and build the source codes by using command:

```
sudo ./install.sh
```

If this way didn't work, try to install bash into your computer. Alternatively, you can build KOMODO manually

```
gfortran -O4 -c src/mod_data.f90
gfortran -O4 -c src/mod_io.f90
gfortran -O4 -c src/mod_xsec.f90
gfortran -O4 -c src/mod_nodal.f90
gfortran -O4 -c src/mod_cmfd.f90
gfortran -O4 -c src/mod_th.f90
gfortran -O4 -c src/mod_trans.f90
gfortran -O4 -c src/mod_control.f90
gfortran -O4 -c src/komodo.f90
gfortran *.o -o komodo
sudo cp KOMODO /usr/bin
```

(NOTE: you need to have admin privilege to run these commands)

These command will create an executable file named `komodo` and copied the executable file to `/usr/bin`. Now, you can run a test using several examples of inputs file in folder [smpl](https://github.com/imronuke/KOMODO/tree/master/smpl) to see if you had built KOMODO properly. You can run KOMODO using command

```
komodo [INPUT_FILE_PATH_NAME]
```

for example, you can run [`IAEA3Ds`](https://github.com/imronuke/KOMODO/blob/master/smpl/static/IAEA3Ds) input by

```
komodo smpl/static/IAEA3Ds
```

If you see `KOMODO EXIT NORMALLY` at the end of terminal output, then congratulations! you have successfully installed KOMODO. By the way, it should take about 0.2 seconds for KOMODO to solve IAEA3Ds problem if you build using gfortran in a typical today's computer. The way to build using Intel fortran is similar, just change `gfortran` with `ifort` in the commands above.

## Building in Mac OS
The way to install KOMODO on Mac OS is similar to installing KOMODO on Ubuntu. First you need to install gfortran to your computer. You can find the information on how to install gfortran for Mac OS from [here](https://gcc.gnu.org/wiki/GFortranBinariesMacOS). And instead of executing this command as you were supposed to do in Ubuntu

```
sudo ./install.sh
```

you shall execute this command in Mac OS

```
sudo ./mac_install.sh
```

The rest of the steps are the same.


## Building in Windows
There are at least two ways installing KOMODO on Windows.

### Using Windows Subsystem for Linux (WSL) in Windows 10
Using this way, you can go to Microsoft Store and install Ubuntu from there for free. Open the Ubuntu app and follow the same steps on installing KOMODO on Ubuntu OS.

### Using min-GW
If you want install quickly on Windows, you can consider using min-GW.

Here is step by step:
1. Download the [KOMODO zip files](https://github.com/imronuke/KOMODO/archive/master.zip) from Github then extract that zip file.
2. Download min-GW from [here](https://github.com/skeeto/w64devkit/releases). Direct link to [64-bit version](https://github.com/skeeto/w64devkit/releases/download/v1.23.0/w64devkit-fortran-1.23.0.zip) and the [32-bit version](https://github.com/skeeto/w64devkit/releases/download/v1.23.0/w64devkit-i686-fortran-1.23.0.zip). Then also unzip min-GW.
3. Run ```w64devkit.exe``` inside min-GW folder, a command prompt will be popped out.
4. In the command prompt, go to the downloaded KOMODO folder. Use ```cd d:``` command to navigate to the D drive.
5. Inside the KOMODO folder, run the following command

    ```
    ./install.sh
    ```
6. An executable ```komodo.exe``` should be generated.
7. Run a test 

```
./komodo.exe smpl\static\IAEA3Ds
```

If you see `KOMODO EXIT NORMALLY` at the end of terminal output, then congratulations! you have successfully installed KOMODO on Windows.
