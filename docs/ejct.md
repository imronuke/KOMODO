---
title: %EJCT
theme: _config.yml
filename: ejct
---

# %EJCT Card

This card us used for transient problems due to control rods insertion or withdrawal.  This card is mandatory if calculation mode is `RODEJECT`.

| %EJCT | Variable | Description | Remarks or examples |
| --- | --- | --- | --- |
| LINE 1 | FBPOS(n) | Final bank n position after ejection and/or insertion (step) | Repeat this line NB times (NB = Number of CR bank) |
|   | TMOVE(n) | When bank n starts move (second) |
|   | BSPEED(n) | Bank n speed (steps/second) |
| LINE 2 | TTOT | Total simulation time (seconds) | Example: `60.  0.1  40.0  1.0` |
|   | TSTEP1 | First time step (seconds) |
|   | TDIV | When to start second time step (seconds) |
|   | TSTEP2 | Second time step (seconds) |
| LINE 3 | IBETA(1:6) | 6-groups delayed neutron fraction | **Not necessary of `%XTAB` card present** |
| LINE 4 | LAMB(1:6) | 6-groups precursor decay constant | **Not necessary of `%XTAB` card present** |
| LINE 5 | VELO(1:NG) | Neutron velocity | **Not necessary of `%XTAB` card present** |

Example:
```
! Rod ejection card
%EJCT
! Final Bank Pos (steps)     Starts Move (seconds)     Speed (steps/seconds)
     228.                     0.0                     2280.0    ! Bank 1     (LINE 1)
     0.                       0.0                     0.0       ! Bank 2
     0.                       0.0                     0.0       ! Bank 3
     228.                     0.0                     0.0       ! Bank 4
     0.                       0.0                     0.0       ! Bank 5
     0.                       0.0                     0.0       ! Bank 6
     0.                       0.0                     0.0       ! Bank 7
5.0 0.005  1.0  0.05            ! (seconds) total time, 1st time step, when to start 2nd time step, 2nd time step
0.0002584 0.00152 0.0013908 0.0030704 0.001102 0.0002584 !beta
0.0128  0.0318  0.119   0.3181  1.4027  3.9286           !decay constants
2.8E7    4.4E5     ! Neutron velocity
```
