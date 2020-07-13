---
title: %PRNT
theme: _config.yml
filename: prnt
---

# %PRNT Card

This card can be used to choose the specific output that users want.

| %PRNT |   |   |   |
| --- | --- | --- | --- |
| LINE 1 | CAPRAD | Radial assembly power distribution print option | 1 = YES<br>0 = NO<br>Default for all = YES |
|   | CAPAXI | Axial assembly power distribution print option |
|   | CAFRAD | Radial Flux Power Distribution |

Example:
```
! PRINT CARD
%PRNT
1 1 0 ! print output
```
