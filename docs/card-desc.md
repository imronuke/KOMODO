---
title: Card Description
theme: _config.yml
filename: card-desc
---

# General Rules

Some general rules for KOMODO inputs:
1.  Input deck is in free-format form with maximum 200 columns
2.	Comments are marked by `!`. Example:
```
! COMPOSITION 1
0.20574  0.00654  0.00415  0.00415  1.0  0.0  0.01462  !Group 1
0.68866  0.04850  0.06099  0.06099  0.0  0.0  0.00000  !Group 2
```

3.	KOMODO input is modular, where it is split into several cards. Cards’ keywords shall be uppercase and marked by `%`. Example:
```
%MODE    ! Calculation mode card
FORWARD
%XSEC    ! Cross section card                                                                                                                                  
2 4      ! Number of groups and number of materials
...
...
%GEOM    ! Geometry card
12 12 2  !nx, ny, nz
...
...
```

4.	Numbers can be repeated using `*` mark. For example
```
10.0 8*20.0  !This line is equivalent to  10.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0
```

5. If you have identical cards for several inputs, it is a good idea to place those cards in separated card files. In the main inputs, you just need to put location of the card files. For example, all cases in NEACRP benchmark have identical cross section data; hence we may place `%XSEC` card in a separated named `neacrp_xsec` which has the following content
```
! XSEC CARD FILE
2  11    ! Number of groups and number of materials
!  sigtr        siga        nu*sigf     kappa*sigf  chi sigs_g1 sigs_g2
5.32058E-02  3.73279E-04  0.00000E+00  0.00000E+00  1.0  0.0  2.64554E-02
3.86406E-01  1.77215E-02  0.00000E+00  0.00000E+00  0.0  0.0  0.00000E+00  !COMP 1
...
...
2.21878E-01  9.71937E-03  6.11444E-03  7.60778E-14  1.0  0.0  1.66175E-02
7.64704E-01  8.85488E-02  1.12635E-01  1.47876E-12  0.0  0.0  0.00000E+00  !COMP 11
```

In the main input, we may refer to the `neacrp_xsec` file by telling KOMODO the location of `neacrp_xsec` preceded by `FILE` keyword
```
! XSEC CARD
%XSEC
FILE /home/imronuke/KOMODO/smpl/static/NEACRP/neacrp_xsec
```


# Card Description

KOMODO has several input cards. Card is a keyword marked with `%`. Each card can be placed arbitrarily in the input deck. Two cards are mandatory for any problems, while the rest are optional and conditional depend on the nature problem being solved. The description of input for each card is explained in this subsection. Following table lists all cards used in KOMODO. You can click each card to get their description


| **No.** | **Cards** | **Description** | **Remark** |
| --- | --- | --- | --- |
| 1. | [`%MODE`](https://imronuke.github.io/KOMODO/mode) | Calculation mode | Mandatory |
| 2. | [`%GEOM`](https://imronuke.github.io/KOMODO/geom) | Geometry of the problem | Mandatory |
| 3. | [`%XSEC`](https://imronuke.github.io/KOMODO/xsec) | Cross Sections | Conditional |
| 4. | [`%CASE`](https://imronuke.github.io/KOMODO/case) | Problem case | Optional |
| 5. | [`%ESRC`](https://imronuke.github.io/KOMODO/esrc) | Extra source | Conditional |
| 6. | [`%ITER`](https://imronuke.github.io/KOMODO/iter) | Iteration Control | Optional |
| 7. | [`%PRNT`](https://imronuke.github.io/KOMODO/prnt) | Output print control | Optional |
| 8. | [`%ADF`](https://imronuke.github.io/KOMODO/adf) | Assembly Discontinuity Factor | Optional |
| 9. | [`%CROD`](https://imronuke.github.io/KOMODO/crod) | Control rods | Conditional |
| 10. | [`%EJCT`](https://imronuke.github.io/KOMODO/ejct) | Control rods ejection and/or insertion | Conditional |
| 11. | [`%FTEM`](https://imronuke.github.io/KOMODO/ftem) | Fuel temperature input card | Conditional |
| 12. | [`%MTEM`](https://imronuke.github.io/KOMODO/mtem) | Moderator/Coolant temperature input card | Conditional |
| 13. | [`%CDEN`](https://imronuke.github.io/KOMODO/cden) | Coolant density input card | Conditional |
| 14. | [`%BCON`](https://imronuke.github.io/KOMODO/bcon) | Boron concentration input card | Conditional |
| 15. | [`%CBCS`](https://imronuke.github.io/KOMODO/cbcs) | Critical boron concentration input card | Conditional |
| 16. | [`%THER`](https://imronuke.github.io/KOMODO/ther) | TH input card | Conditional |
| 17. | [`%XTAB`](https://imronuke.github.io/KOMODO/xtab) | XSEC library for branch calculations | Conditional |
| 18. | [`%KERN`](https://imronuke.github.io/KOMODO/kern) | Nodal kernel options | Optional |
| 19. | [`%EXTR`](https://imronuke.github.io/KOMODO/extr) | Exponential flux transformation option card for transient problem | Optional |
| 20. | [`%THET`](https://imronuke.github.io/KOMODO/thet) | Used to set theta value for transient problem | Optional |
