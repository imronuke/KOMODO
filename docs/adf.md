---
title: %ADF
theme: _config.yml
filename: adf
---

# %ADF Card

ADF card can be incorporated into ADPRES input, if any, to make the solution more accurate.

| %ADF | Variable | Description | Remarks or examples |
| --- | --- | --- | --- |
| LINE 1 | DC(1) | East side discontinuity factor | Repeat LINE 2 NG times. And again repeat this input segment NMAT times.(See example below) |
|   | DC(2) | West side discontinuity factor |
|   | DC(3) | North side discontinuity factor |
|   | DC(4) | South side discontinuity factor |
|   | DC(5) | Bottom side discontinuity factor |
|   | DC(6) | Top side discontinuity factor |
| LINE 2 | ROT | Direction for ADF rotation | This line is optional. Necessary if the value of discontinuity factor is not fully symmetric. Repeat this line followed by LINE 3 as many as desired until zero number (LINE 5) is entered.<br>1 = 90 degree counter clockwise<br>2 = 180 degree counter clock wise<br>3 = 270 degree counter clockwise |
| LINE 3 | X1 | Start assembly position in X-direction | This line follows LINE 2 which tells the position of assembly being rotated. Repeat this line as many as desired until zero numbers (LINE 4). Followings are value limits for these line<br>1 \&lt;= X1 \&lt;= NX;<br>1 \&lt;= X2 \&lt;= NX;<br>1 \&lt;= Y1 \&lt;= NY;<br>1 \&lt;= Y2 \&lt;= NY;<br>1 \&lt;= Z1 \&lt;= NZ;<br>1 \&lt;= Z2 \&lt;= NZ;<br>X1 \&lt;= X2; Y1 \&lt;= Y2; Z1 \&lt;= Z2 |
|   | X2 | End assembly position in X-direction |
|   | Y1 | Start assembly position in Y-direction |
|   | Y2 | End assembly position in Y-direction |
|   | Z1 | Start assembly position in Z-direction |
|   | Z2 | End assembly position in Z-direction |
| LINE 4 | 0 | Zero numbers entered to end X1 through Z2 |  |
|   | 0 |
|   | 0 |
|   | 0 |
|   | 0 |
|   | 0 |
| LINE 5 | 0 | Zero number to end ROT |
| LINE 6 | ZP | ADF print option | This line is optional. |

Example:
```
! Assembly Discontinuity Factors Card
%ADF
!g1 -> east  west  north  south
!g2 -> east  west  north  south
! COMPOSITION 1
  0.9966 0.9288 0.9966  0.9288  1.0000  1.0000   ! LINE 1
  1.1332 1.6570 1.1332  1.6570  1.0000  1.0000   ! LINE 1
! COMPOSITION 2
  1.0787 0.8423 1.0787  0.8423  1.0000  1.0000   ! LINE 1
  1.6423 0.6809 1.6423  0.6809  1.0000  1.0000   ! LINE 1
! COMPOSITION 3
  0.9989 0.9114 0.9989  0.9114  1.0000  1.0000   ! LINE 1
  1.1664 1.5805 1.1664  1.5805  1.0000  1.0000   ! LINE 1
! COMPOSITION 4
1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   ! LINE 1
1.0000  1.0000  1.0000  1.0000  1.0000  1.0000   ! LINE 1
! ADF ROTATION DUE TO DIAGONALLY SYMMETRIC ADFs
 1     ! 90 degree counter clock-wise            ! LINE 2
 2  2  1  1  1  2                                ! LINE 3 (and next rows)
 2  2  3  3  1  2
 2  2  5  5  1  2
 2  2  7  7  1  2
 2  2  9  9  1  2
 2  2 11 11  1  2
 4  4  1  1  1  2
 4  4  3  3  1  2
 4  4  5  5  1  2
 4  4  7  7  1  2
 4  4  9  9  1  2
 6  6  1  1  1  2
 6  6  3  3  1  2
 6  6  5  5  1  2
 6  6  7  7  1  2
 6  6  9  9  1  2
 8  8  1  1  1  2
 8  8  3  3  1  2
 8  8  5  5  1  2
 8  8  7  7  1  2
10 10  1  1  1  2
10 10  3  3  1  2
 0  0  0  0  0  0                                ! LINE 4
 2      ! 180 degree counter clock-wise
 2  2  2  2  1  2                                ! LINE 3 (and next rows)
 2  2  4  4  1  2
 2  2  6  6  1  2
 2  2  8  8  1  2
 2  2 10 10  1  2
 4  4  2  2  1  2
 4  4  4  4  1  2
 4  4  6  6  1  2
 4  4  8  8  1  2
 4  4 10 10  1  2
 6  6  2  2  1  2
 6  6  4  4  1  2
 6  6  6  6  1  2
 6  6  8  8  1  2
 8  8  2  2  1  2
 8  8  4  4  1  2
 8  8  6  6  1  2
 8  8  8  8  1  2
10 10  2  2  1  2
10 10  4  4  1  2
 0  0  0  0  0  0                                ! LINE 4
 3     ! 270 degree counter clock-wise
 1  1  2  2  1  2                                ! LINE 3 (and next rows)
 1  1  4  4  1  2
 1  1  6  6  1  2
 1  1  8  8  1  2
 1  1 10 10  1  2
 3  3  2  2  1  2
 3  3  4  4  1  2
 3  3  6  6  1  2
 3  3  8  8  1  2
 3  3 10 10  1  2
 5  5  2  2  1  2
 5  5  4  4  1  2
 5  5  6  6  1  2
 5  5  8  8  1  2
 7  7  2  2  1  2
 7  7  4  4  1  2
 7  7  6  6  1  2
 7  7  8  8  1  2
 9  9  2  2  1  2
 9  9  4  4  1  2
 9  9  6  6  1  2
11 11  2  2  1  2
 0  0  0  0  0  0                                ! LINE 4
 0   ! ADF INPUTS END                            ! LINE 5
 1   ! Print option                              ! LINE 6
```
