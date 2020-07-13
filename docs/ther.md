---
title: %THER
theme: _config.yml
filename: ther
---

# %THER Card

This card is used to set the T-H parameters. It also activates T-H feedback. When  `%THER` card present, either conditions must be satisfied
1. `%FTEM` and one or more of these cards shall also present: `%MTEM` and `%CDEN`
2. `%XTAB` shall be present

| `%THER` | Variable | Description | Remarks |
| --- | --- | --- | --- |
| LINE 1 | PPOW | Percent reactor power | 0.0 < PPOW <= 100.0 |
| LINE 2 | POW | Reactor thermal power for given geometry in `%GEOM` card in Watt |  |
| LINE 3 | TIN | Coolant inlet temperature in Kelvin |  |
|        | CMFLOW | Coolant mass flow rate in kg/s    |  |
| LINE 4 | RF | Fuel meat (UO2) radius in meter |  |
|        | TG | Gap thickness in meter    |  |
|        | TC | Cladding thickness in meter    |  |
|        | PPITCH | Fuel pin pitch in meter    |  |
| LINE 5 | NFPIN | Number of fuel pin for each assembly |  |
|        | NGT | Number of guide tubes    |  |
| LINE 6 | CF | Fraction of heat deposited in the coolant |  |

Example:
```
! THERMAL-HYDRAULIC CARD
%THER
100.                                ! Percent power in %
891.25e6                            ! Reactor thermal power for quarter geometry in Watt
560.  82.121243523                  ! Inlet coolant temp. (Kelvin) and Fuel Assembly Mass flow rate (kg/s)
3.951E-03  5.9E-05  5.73E-04  1.26E-2 ! Fuel meat rad., gap thinkness, cladding thinkness and pin picth (m)
264  25                               ! Number of fuel pin and guide tubes
0.0                                   ! FRACTION OF HEAT DEPOSITED IN COOLANT
1
```
