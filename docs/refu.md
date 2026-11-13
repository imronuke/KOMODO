---
title: %REFU
theme: _config.yml
filename: refu
---

# %REFU Card

`%REFU` loads assembly-wise node state from a previous `%CYCL` cycle file and assigns it to the current core map. `%DLIB` is required. The card can be used with `%DEPL` to continue depletion or with a regular DLIB calculation.

The first line names the prior-cycle file. It is followed by `NY` assembly-map rows, each containing exactly `NX` entries. Rows are written from `+Y` to `-Y`, with columns from `-X` to `+X`:

```text
%REFU
cycle1.cyl
0 0 0 0 0
0 A01 XXX B01 0
0 C01 D01 E01 0
0 0 0 0 0
```

Each nonzero ID other than `XXX` identifies one unique assembly ID in the prior-cycle file. KOMODO compares every node material in that source assembly with the material at the corresponding destination node in `%GEOM`. If they match, KOMODO transfers that node's burnup and history boron, moderator temperature, and fuel temperature. The same node layout and orientation are required; assembly rotation is not supported.

`XXX` marks fresh fuel. KOMODO skips source-assembly material matching and keeps the initial burnup and reference history conditions specified by `%DLIB` for the destination materials. `0` denotes a nonfuel or outside-core position and must match the current `%GEOM` fuel layout.

The prior-cycle geometry dimensions, axial node widths, and fuel assembly-local node widths must match the current geometry. A prior-cycle assembly ID cannot be used more than once. Assemblies absent from the new map are treated as discharged.
