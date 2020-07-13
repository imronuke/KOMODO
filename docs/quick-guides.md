---
title: Quick Guides
theme: _config.yml
filename: quick-guides
---

# Quick Guides
## Writing Input
ADPRES input is designed to be self-explanatory. It has several input cards marked by `%`, for example: `%mode`, `%geom`, `%xsec`, and so on. Some cards are mandatory for any problems. While some cards are conditional depending on the problems at hand and some cards are optional. Comments are marked by `!`. For example, the following is the [IAEA3D input](https://github.com/imronuke/ADPRES/blob/master/smpl/static/IAEA3Ds), where you can find its specification [here](https://engineering.purdue.edu/PARCS/Code/TestSuite/CalculationMode/StandAloneMode/Eigenvalue/IAEA3DPWR).

```
! IAEA3D input data
! NODE SIZE = 10 cm
! PARCS K-EFF  : 1.029096
! ADPRES K-EFF : 1.029082 (ERROR = 1.4 PCM)

! Mode card
%MODE
FORWARD

! Case card
%CASE
IAEA3D
10 CM NODE SIZE

! Cross-sections card
%XSEC
2  5    ! Number of groups and number of materials
! sigtr   siga   nu*sigf sigf   chi   sigs_g1  sigs_g2
0.222222  0.010  0.000  0.000    1.0   0.1922   0.020
0.833333  0.080  0.135  0.135    0.0   0.000    0.7533   ! MAT1 : Outer Fuel
0.222222  0.010  0.000  0.000    1.0   0.1922   0.020
0.833333  0.085  0.135  0.135    0.0   0.000    0.7483   ! MAT2 : Inner Fuel
0.222222  0.0100 0.000  0.000    1.0   0.1922   0.020
0.833333  0.1300 0.135  0.135    0.0   0.000    0.7033   ! MAT3 : Inner Fuel + Control Rod
0.166667  0.000  0.000  0.000    0.0   0.1267   0.040
1.111111  0.010  0.000  0.000    0.0   0.000    1.1011   ! MAT4 : Reflector
0.166667  0.000  0.000  0.000    0.0   0.000    0.040
1.111111  0.055  0.000  0.000    0.0   0.000    0.000    ! MAT5 : Reflector + Control Rod
%GEOM
9 9 19         ! number of assembly in x, y, z directions
10.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0    !x-direction assembly size in cm
1      8    8    8    8    8    8    8    8     !x-direction assembly divided into 2 (10 cm each)
20.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0 10.0    !y-direction assembly size in cm
8      8    8    8    8    8    8    8    1     !y-direction assembly divided into 2 (10 cm each)
19*20.0                                         !z-direction assembly  in cm
19*1                                            !z-direction nodal is not divided
4                                               !np number of planar type
1  13*2  4*3  4                                 !planar assignment (from bottom to top)
! Planar_type_1 (Bottom Reflector)
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  0
  4  4  4  4  4  4  4  4  0
  4  4  4  4  4  4  4  0  0
  4  4  4  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Planar_type_2 (Fuel)
  3  2  2  2  3  2  2  1  4
  2  2  2  2  2  2  2  1  4
  2  2  2  2  2  2  1  1  4
  2  2  2  2  2  2  1  4  4
  3  2  2  2  3  1  1  4  0
  2  2  2  2  1  1  4  4  0
  2  2  1  1  1  4  4  0  0
  1  1  1  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Planar_type_3 (Fuel+Partial Control Rods)
  3  2  2  2  3  2  2  1  4
  2  2  2  2  2  2  2  1  4
  2  2  3  2  2  2  1  1  4
  2  2  2  2  2  2  1  4  4
  3  2  2  2  3  1  1  4  0
  2  2  2  2  1  1  4  4  0
  2  2  1  1  1  4  4  0  0
  1  1  1  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Planar_type_4 (Top reflectors)
  5  4  4  4  5  4  4  4  4
  4  4  4  4  4  4  4  4  4
  4  4  5  4  4  4  4  4  4
  4  4  4  4  4  4  4  4  4
  5  4  4  4  5  4  4  4  0
  4  4  4  4  4  4  4  4  0
  4  4  4  4  4  4  4  0  0
  4  4  4  4  4  4  0  0  0
  4  4  4  4  0  0  0  0  0
! Boundary conditions
! 0 = zero-flux
! 1 = zero-incoming current
! 2 = reflective
(east), (west), (north), (south), (bottom), (top)
   1       2       2        1        1        1

! NOTE: Writing 19*20.0 is equivalent to write 20.0 nineteen times in a row
```

In the above example, there are
1. Two mandatory cards : `%MODE` and `%GEOM`
2. One conditional card: `%XSEC`
3. One optional card   : `%CASE`

You find the detailed description for each card [here](https://imronuke.github.io/ADPRES/card-desc), but we will explain them briefly now
### [%MODE Card](https://imronuke.github.io/ADPRES/mode)
This is the mode of ADPRES calculation. Since here we want to calculate static forward calculation (eigenvalue problem) the calculation mode is `FORWARD`. The detailed description of this card is [here](https://imronuke.github.io/ADPRES/mode).
### [%CASE Card](https://imronuke.github.io/ADPRES/case)
This card is optional. This describes the problem at hand. The detailed description of this card is [here](https://imronuke.github.io/ADPRES/case).
### [%XSEC Card](https://imronuke.github.io/ADPRES/xsec)
This card is conditional, necessary only if `%XTAB` card is not present. This card tells ADPRES the cross sections data for the problem. The cross section data must be given for each group and for each material as shown in the example. The description of the cross sections data can be seen in the comments. The detailed description of this card is [here](https://imronuke.github.io/ADPRES/xsec).
### [%GEOM Card](https://imronuke.github.io/ADPRES/geom)
This card is describes the geometry of the problem. It quite similar to other reactor core simulator which you can easily understand if you have background on nuclear engineering. The description of the inputs given in the comments. The detailed description of this card is [here](https://imronuke.github.io/ADPRES/geom).


## Running a Test
In Linux or other Unix based OS, you can run ADPRES using command

```
adpres [INPUT_FILE_PATH_NAME]
```

for example, you can run [`IAEA3Ds`](https://github.com/imronuke/ADPRES/blob/master/smpl/static/IAEA3Ds) input by

```
adpres /home/imronuke/smpl/static/IAEA3Ds
```

While in Windows, for example, you can run as follow

```
adpres C:\Users\imronuke\Downloads\ADPRES-master\smpl\static\IAEA3Ds
```


## Reading Output

After you run a test, you should see in the summary of the output in terminal as follow

```
           ###########################################################
           #                     ADPRES 1.2                          #
           #        ABU DHABI POLYTECHNIC REACTOR SIMULATOR          #
           ###########################################################


  CALCULATION MODE : FORWARD CALCULATION                                         

  CASE ID : IAEA3D                                                                                              
  10 CM NODE SIZE                                                                                     

  NODAL KERNEL  : SEMI-ANALYTIC NODAL METHOD

  reading input ... done


  ==============================================================================
                           CALCULATION RESULTS
  ==============================================================================

  Itr     k-eff     Fis.Src Error   Inner Error
 ----------------------------------------------------
    1     0.981424    5.47871E-01    8.55259E+03
    2     1.001319    2.41976E-01    8.56379E+00
    3     1.009804    1.67297E-01    6.27992E-01
    4     1.013500    1.33833E-01    1.44559E-01
     ...FISSION SOURCE EXTRAPOLATED...
    5     1.026980    2.65156E+00    1.31348E-01
                       .
                       .
                       .
   15     1.028565    5.05638E-01    6.73943E-02
   16     1.028253    6.31685E-02    4.20531E-01
   17     1.028196    3.63469E-02    4.93252E-01
   18     1.028082    2.21127E-02    2.48688E-01
   19     1.028020    1.29756E-02    6.24972E-02
     ...FISSION SOURCE EXTRAPOLATED...
   20     1.027400    3.54750E-01    2.34432E-02
   21     1.027655    6.04390E-02    3.09183E-01
     .....NODAL COUPLING UPDATED.....
MAX. CHANGE IN NODAL COUPLING COEF.=  3.16843E-01 AT NODE I =  6, J =  4, K = 19
   22     1.028084    5.88420E-02    1.59061E-01
   23     1.032373    1.30872E-01    4.33336E-01
                      .
                      .
                      .
  127     1.029082    1.20168E-05    8.28493E-05
  128     1.029082    9.48481E-06    2.04727E-05
  129     1.029082    7.46358E-06    9.66768E-06

  MULTIPLICATION EFFECTIVE (K-EFF) =  1.029082


  CPU time breakdown in seconds
    Input reading time   :    0.0078  ( 5.0%)
    XSEC processing time :    0.0003  ( 0.2%)
    CMFD time            :    0.0893  (57.1%)
    Nodal update time    :    0.0590  (37.7%)
    T-H time             :    0.0000  ( 0.0%)
    ------------------------------------------
    Total time           :    0.1565

  ADPRES EXIT NORMALLY
  ```

  If you get `ADPRES EXIT NORMALLY` in the end of the terminal output, it means you successfully run ADPRES. Since it is a forward (eigenvalue) problem, you will see the outer iterations as they evolve and you will see also the effective multiplication factor as well as CPU time breakdown. The detailed output, such as radial and axial power distribution, can be found in the same file name as input but with an `.out` extension.

  It is always a good idea to see this output file to ensure that you had written input correctly. ADPRES echoes your input and redescribe the input to make sure that this is the problem you want to solve.  ADPRES may run well without any error but it gives you wrong results. This might happen if the input is not consistent with the problem specification.
