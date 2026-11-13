---
title: %DEPL
theme: _config.yml
filename: depl
---

# %DEPL Card

This card gives the operating schedule for depletion intervals. Both `%DLIB` and `%THER` must be present. The card is valid only in `FORWARD` and `BCSEARCH` modes. In `BCSEARCH` mode, `%CBCS` is required.

| `%DEPL` | Variable | Description |
| --- | --- | --- |
| LINE 1 | N-INTERVALS | Number of operating intervals |
| NEXT N-INTERVALS LINES | DAYS | Interval duration in days, followed by reactor power level in percent and one control-bank step position for each bank defined in `%CROD` |

Example:

```text
%DEPL
2
100.0 100.0 120 120
150.0 80.0 100 100
```

Each schedule row has the form `DAYS POWER-LEVEL-% BANK-1-STEP ... BANK-NB-STEP`. If `%CROD` is absent, each row contains only duration and power level. The schedule power level overrides the percent power on line 1 of `%THER`; line 2 of `%THER` remains the rated full thermal power used to calculate interval MWth and burnup. The power level must be greater than 0 and no more than 100 percent.

Thermal feedback supplies fuel and moderator temperatures and moderator density throughout the schedule.
