---
title: %ITER
theme: _config.yml
filename: iter
---

# %ITER Card

This card can be used to control the iterations in ADPRES calculation

| %ITER | Variable | Description | Remarks |
| --- | --- | --- | --- |
| LINE 1 | NOUT  | Maximum number of outer iteration | Default = 500 |
|        | NIN   | Maximum number of inner iteration per an outer iteration | Default = 2 |
|        | SERC  | Fission source error criteria (relative error) | Default = 1.e-5 |
|        | FERC  | Flux error criteria (relative error) | Default = 1.e-5 |
|        | NAC   | Outer iteration fission source extrapolation interval | Default = 5 |
|        | NUPD  | Nodal update interval through outer iteration | Default = (NX+NY+NZ)/2.5. Effective only for ANM and PNM [nodal kernel](https://imronuke.github.io/ADPRES/kern) |
|        | TH_NITER  | Maximum number of T-H iteration | Default = 30. Effective only if [`%THER`](https://imronuke.github.io/ADPRES/ther) present  |
|        | NTH  | Maximum number of outer iteration per T-H iteration | Default = 20. Effective only if [`%THER`](https://imronuke.github.io/ADPRES/ther) present |

Note: TH iteration is non-linear iteration to determine TH parameters

Example:
```
! Iteration control card
%ITER
! Set maximum number of outer iteration = 1200
! Set maximum number of inner iteration = 3
! Set fission source error criteria = 1.e-5
! Set flux error criteria = 1.e-5
! Set fission source extrapolation interval every 25 outer iterations
! Set nodal update interval every 40 outer iterations
! Set maximum number of TH iterations = 20
! Set maximum number outer iterations for each TH iteration is 80
1200 3 1.e-5 1.e-5 25 40 20 80
```
