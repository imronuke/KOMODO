---
title: %XTAB
theme: _config.yml
filename: xtab
---

# %XTAB Card

This card is used for branch calculations. The format library used in ADPRES is in tabular form similar to the ones used in [UO2/MOX Transient Benchmark](https://www.oecd-nea.org/science/wprs/MOX-UOX-transients/benchmark_documents/final_report/mox-uo2-bench.pdf). An example of this library can be seen [here](https://github.com/imronuke/ADPRES/blob/master/smpl/xsec/SERPENT_CMM/m40.tab). This card is conditional : either `%XSEC` or `XTAB` card must be present

| `%XTAB` | Variable | Description | Remarks |
| --- | --- | --- | --- |
| LINE 1 | NG | Number of groups |  |
|        | NMAT | Number materials |
| LINE 2 | FNAME | The path of the library | Repeat  this line NMAT times.(See example below) |
|  | FNUM | Number of the material inside the library FNAME |

Example:
```
! XTAB CARD
%XTAB
2  19    ! Number of groups and number of materials
! XTAB files for each material
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u42.tab      1  !MAT  1 : u42u  0.15 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u42.tab      2  !MAT  2 : u42u  17.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u42.tab      4  !MAT  3 : u42u  22.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u42.tab      5  !MAT  4 : u42u  32.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u42.tab      6  !MAT  5 : u42u  35.0 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u42.tab      7  !MAT  6 : u42u  37.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u45.tab      1  !MAT  7 : u45u  0.15 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u45.tab      2  !MAT  8 : u45u  17.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u45.tab      3  !MAT  9 : u45u  20.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u45.tab      5  !MAT  10: u45u  32.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/u45.tab      7  !MAT  11: u45u  37.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/m40.tab      1  !MAT  12: m40m  0.15 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/m40.tab      4  !MAT  13: m40m  22.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/m40.tab      7  !MAT  14: m40m  37.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/m43.tab      1  !MAT  15: m43m  0.15 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/m43.tab      2  !MAT  16: m43m  17.5 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/m43.tab      6  !MAT  17: m43m  35.0 GWD/TON
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/refr.tab     1  !MAT  18: RADIAL REFLECTOR
/home/imronuke/ADPRES/smpl/xsec/SERPENT_CMM/refa.tab     1  !MAT  19: AXIAL REFLECTOR
1
```
