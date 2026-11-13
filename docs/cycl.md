---
title: %CYCL
theme: _config.yml
filename: cycl
---

# %CYCL Card

`%CYCL` saves the completed depletion state for later fuel shuffling. It is valid only with `%DEPL` and `%DLIB`. Use the saved cycle file with [`%REFU`](https://imronuke.github.io/KOMODO/refu) to assign prior-cycle assemblies to a new core map.

The card contains an assembly map followed by the output filename:

```text
%CYCL
! One row per assembly row, from +Y to -Y; columns are from -X to +X
0 0 A01 B01 0
0 C01 D01 E01 0
0 0 0 0 0
cycle1.cyl
```

Each map row must contain exactly `NX` entries, and the map must have `NY` rows, matching `%GEOM`. Every nonzero assembly ID must be unique. `0` may repeat and denotes a non-fuel or outside-core position. Every position containing fuel must have a nonzero ID, and every nonzero ID must identify a position containing fuel.

Fuel assemblies must have the same node subdivisions and node widths in the x and y directions. Their local orientation must match the `%GEOM` orientation; rotating assemblies is not supported. The global axial node layout is shared by the geometry.

After all depletion intervals complete, KOMODO prints the assembly IDs and the axially averaged burnup map to the `.out` log and writes the cycle state to the specified file. The map burnup is HM-mass weighted and reported in `MWd/kgHM`. The cycle file stores the common assembly-local x/y node widths and global axial node widths, plus each node's assembly and local node indices, axial index, DLIB material name, burnup, history boron, history moderator temperature, and history fuel temperature. Non-fuel node burnup is stored as zero. A cycle file is written only when the full schedule completes.
