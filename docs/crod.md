---
title: %crod
theme: _config.yml
filename: crod
---

# %CROD Card

If the problems have control rods inserted, users can use this card. This card is mandatory for `RODEJECT` calculation mode.

| %CROD | Variable | Description | Remarks or examples |
| --- | --- | --- | --- |
| LINE 1 | NB | Number of CR banks |   |
|   | NSTEP | Maximum number of steps |
| LINE 2 | POS0 | Zero step position (cm from bottom) |   |
|   | SSIZE | step size (cm/step) |
| LINE 3 | BPOS(1:NB) | Control Rod Bank position (step) **0 step means full inserted** |  |
| LINE 4 | BMAP(1:NX) | Control Rod Bank Map | Repeat this line NY times (see example the input below) |
| LINE 5 | DISGTR(g) | Macroscopic Cross Section changes due to control rods insertion | Repeat LINE 2 NG times. And again repeat this input segment NMAT times. **This line is not necessary if `%XTAB` card present** |
|   | DSIGA(g) |
|   | DNUF(g) |
|   | DSIGF(g) |
|   | CHI(g) |
|   | DSIGS(g,1:NG) |

Example:
```
! CONTROL CARD
%CROD
7    228                            ! Number of CR banks and max number of banks
37.7 1.5942237                      ! Zero step pos. (cm) and cm/step (total 228 steps)
0.   0.  0.  228.   0.   0.  0.     ! CR Bank pos. (0=fully inserted, 228=fully withdrawn)
 1  0  2  0  0  0  3  0  0          ! (LINE 4)
 0  4  0  0  0  6  0  0  0
 2  0  5  0  6  0  6  0  0
 0  0  0  4  0  0  0  0  0
 0  0  6  0  7  0  0  0  0
 0  6  0  0  0  0  0  0  0
 3  0  6  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
 0  0  0  0  0  0  0  0  0
! CX changes
!  sigtr        siga        nu*sigf        kappa*sigf      sigs_g1   sigs_g2
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 1
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 2
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 3
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 4
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 5
 3.74092E-03  2.42926E-03  -1.22634E-04  -1.47557E-15  0.0  -3.14239E-03
-1.67503E-02  2.56478E-02  -3.28086E-03  -4.30444E-14  0.0   0.00000E+00    !COMP 6
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 7
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 8
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 9
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 10
 3.73220E-03  2.47770E-03  -1.02786E-04  -1.21448E-15  0.0  -3.19253E-03
-2.19926E-02  2.55875E-02  -2.82319E-03  -3.70238E-14  0.0   0.00000E+00    !COMP 11
```
