---
title: %DLIB
theme: _config.yml
filename: dlib
---

# %DLIB Card

This card selects a depletion library and maps its compositions to KOMODO material IDs. The energy-group count is read from the library's `GROUPS` entry. Exactly one of `%XSEC`, `%XTAB`, or `%DLIB` must be present as the cross-section source. If `%DEPL` is present, `%DLIB` is required.

| `%DLIB` | Variable | Description |
| --- | --- | --- |
| LINE 1 | FNAME | Path to the `.dlib` library |
| LINE 2 | NMAT | Number of material IDs to map |
| NEXT NMAT LINES | COMPOSITION-NAME | Library material name, followed by the initial burnup in `MWd/kgHM`; rows map to material IDs in order, starting with MAT-ID 1 |

Example:

```text
%DLIB
casmo/apr1400.dlib
3
br 0.0
tr 0.0
a0 12.5
```

Names are matched case-insensitively against the `MATERIAL` names in the library. Each name must be unique in the library and present in its `MATERIAL` entries; numeric composition IDs are not accepted in these rows.

For each node, the history boron concentration and fuel and moderator temperatures start at the material's reference conditions. This assumes the specified initial burnup accumulated under those reference conditions. After each interval, each history value is updated using the node's burnup increment:

```text
H_new = (BU_old * H_old + delta_BU * current_condition) / BU_new
```

The updated history values are used for the next cross-section lookup. Nodes without a burnup increment retain their previous history values.

DLIB lookups linearly extrapolate current and history temperature and boron conditions using the nearest two stored points, including the condition axes in reflector interaction data. If a condition dimension has only one stored point, its correction is held constant outside that point. Burnup remains strictly bounded; requests outside a curve's burnup range are errors. Floating-point values within `1e-8` of a stored endpoint are treated as that endpoint.

ROD curves are applied to fuel nodes only. If a material has one ROD curve, KOMODO uses it for every control-rod bank; if it has multiple ROD curves, the curve label must match the bank label (for example, `CR3`).
