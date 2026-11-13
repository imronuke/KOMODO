---
title: %HDF5
theme: _config.yml
filename: hdf5
---

# %HDF5 Card

This card writes the reactor state to an HDF5 file for post-processing, for example with Python (`h5py`). KOMODO must be built with the CMake option `HDF5_OUTPUT=ON` (the default); otherwise the card stops the run with an error.

| `%HDF5` | Variable | Description | Remark |
| --- | --- | --- | --- |
| LINE 1 | H5FILE | Name of the HDF5 output file | Optional. If omitted, the file is named after the input file with the `.h5` extension |

Example with the default file name (input file `depletion` gives `depletion.h5`):
```
%HDF5
```

Example with a user-defined file name:
```
%HDF5
apr1400_cycle1.h5
```

## File Layout

The file contains three top-level groups:

| Group | Content |
| --- | --- |
| `/metadata` | Schema version, input file name, calculation mode, number of groups and nodes, `nx`, `ny`, `nz`, `nxx`, `nyy`, `nzz`, and flags for T-H and transient data. With `%THER`, also the full, percent, and operating power |
| `/geometry` | Node indices `ix`, `iy`, `iz`, node `material`, `node_volume_cm3`, and node widths `x_width_cm`, `y_width_cm`, `z_width_cm`. With `%THER`, also the fuel radial mesh positions |
| `/steps/NNNNNN` | One group per saved state, numbered with six digits starting from `000000` |

A state is saved once for static calculations (`FORWARD`, `ADJOINT`, `FIXEDSRC`, `BCSEARCH`), at the initial state and every time step for transients, and at the start of every depletion interval for `%DEPL`.

Each step group has a `global` and a `node` subgroup:

| Subgroup | Datasets |
| --- | --- |
| `global` | `time_s`, `normalization`, units of flux, reaction rate, and power density, `reactivity` and `reactivity_pcm` (not for `FIXEDSRC`; plus `reactivity_dollar` for transients). Depending on the problem, also `power_W`, `boron_concentration_ppm`, `control_rod_bank_position_step` and `control_rod_bank_position_cm`, and maximum and average fuel temperature, moderator temperature, and coolant density |
| `node` | `flux` (node x group), `nu_fission_rate` (node x group), and `power_density`. With `%THER`, also `fuel_temperature_K`, `moderator_temperature_K`, `coolant_density_g_cm3`, `fuel_centerline_temperature_K`, and `clad_surface_temperature_K`. With `%DEPL`, also `burnup_MWd_kgHM` |

With `%THER`, flux and power density are absolute (n/cm<sup>2</sup>/s and kW/cm<sup>3</sup>), normalized to the operating power. Without `%THER`, they are relative.
