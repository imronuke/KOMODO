---
title: Card Description
theme: _config.yml
filename: card-desc
---

# General Rules

Some general rules for ADPRES inputs:
1.  Input deck is in free-format form with maximum 200 columns
2.	Comments are marked by `!`. Example:
```
! COMPOSITION 1
0.20574  0.00654  0.00415  0.00415  1.0  0.0  0.01462  !Group 1
0.68866  0.04850  0.06099  0.06099  0.0  0.0  0.00000  !Group 2
```

3.	ADPRES input is modular, where it is broken into several cards. Cards’ keywords shall be uppercase and marked by `%`. Example:
```
%MODE
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
10.0 8*20.0  !is equivalent to  10.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0 20.0
```


# Card Description

ADPRES has several input cards. Card is a keyword marked with `%`. Each card can be placed arbitrarily in the input deck. Two cards are mandatory for any problems, while the rest are optional and conditional depend on the nature problem being solved. The description of input for each card is explained in this subsection. Following table lists all cards used in ADPRES. You can click each cards to get their description


| **No.** | **Cards** | **Description** | **Remark** |
| --- | --- | --- | --- |
| 1. | [`%MODE`](https://imronuke.github.io/ADPRES/mode) | Calculation mode | Mandatory |
| 2. | [`%GEOM`](https://imronuke.github.io/ADPRES/geom) | Geometry of the problem | Mandatory |
| 3. | [`%XSEC`](https://imronuke.github.io/ADPRES/xsec) | Cross Sections | Conditional |
| 4. | [`%CASE`](https://imronuke.github.io/ADPRES/case) | Problem case | Optional |
| 5. | [`%ESRC`](https://imronuke.github.io/ADPRES/esrc) | Extra source | Conditional |
| 7. | [`%ITER`](https://imronuke.github.io/ADPRES/iter) | Iteration Control | Optional |
| 8. | [`%PRNT`](https://imronuke.github.io/ADPRES/prnt) | Output print control | Optional |
| 9. | [`%ADF`](https://imronuke.github.io/ADPRES/adf) | Assembly Discontinuity Factor | Optional |
| 10. | [`%CROD`](https://imronuke.github.io/ADPRES/crod) | Control rods | Conditional |
| 11. | [`%EJCT`](https://imronuke.github.io/ADPRES/ejct) | Control rods ejection and/or insertion | Conditional |
| 12. | [`%FTEM`](https://imronuke.github.io/ADPRES/ftem) | Fuel temperature input card | Conditional |
| 13. | [`%MTEM`](https://imronuke.github.io/ADPRES/mtem) | Moderator/Coolant temperature input card | Conditional |
| 14. | [`%CDEN`](https://imronuke.github.io/ADPRES/cden) | Coolant density input card | Conditional |
| 15. | [`%BCON`](https://imronuke.github.io/ADPRES/bcon) | Boron concentration input card | Conditional |
| 16. | [`%CBCS`](https://imronuke.github.io/ADPRES/cbcs) | Critical boron concentration input card | Conditional |
| 17. | [`%THER`](https://imronuke.github.io/ADPRES/ther) | TH input card | Conditional |
| 18. | [`%XTAB`](https://imronuke.github.io/ADPRES/xtab) | XSEC library for branch calculations | Conditional |
| 19. | [`%KERN`](https://imronuke.github.io/ADPRES/kern) | Nodal kernel options | Optional |
| 20. | [`%EXTR`](https://imronuke.github.io/ADPRES/extr) | Exponential flux transformation option card for transient problem | Optional |
| 21. | [`%THET`](https://imronuke.github.io/ADPRES/thet) | Used to set theta value for transient problem | Optional |
