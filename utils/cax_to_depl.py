'''
To convert CASMO CAX files to KOMODO depletion library (dlib).
Currently only work for CASMO input with S3C 'HOT'
For reflectors, at least two temperature points must be available. Try use:

 DEP  0                   * CARD INSERTED BY S3C PREPROCESSOR
 WRI 0/'RES' NLI
 WRI 0/'OUT'
 STA  * (Card added by SCANIN)
 * S3C 'HOT'
 * PWR - RADIAL REFLECTOR DATA -
 COE ,/0001, 0.0
 ROD='OUT' BOR=0.1 1000.0 2000.0 TMO= 560.4 600.4 GAM=1000
 + REF 0.27940 2.2225 18.2753 / 'MOD' 'CRS' 'MOD' /// 'RAD'

'''

from __future__ import annotations

import argparse
from collections.abc import Sequence
from dataclasses import dataclass, replace
import math
from pathlib import Path
import re
import shlex

if __package__:
    from .cax_to_tab import (
        CaxLibrary,
        ScalarTable,
        State,
        branch_delta,
        build_library,
        classify_state,
        has_reflector_record,
        is_direct_rod_state,
        parse_condition_line,
        parse_record_rows,
        read_cax_states,
    )
else:
    from cax_to_tab import (
        CaxLibrary,
        ScalarTable,
        State,
        branch_delta,
        build_library,
        classify_state,
        has_reflector_record,
        is_direct_rod_state,
        parse_condition_line,
        parse_record_rows,
        read_cax_states,
    )


_CURRENT_BRANCHES = (
    ("BORON", "current_bor", "base_bor"),
    ("MODERATOR_TEMPERATURE", "current_tmo", "base_tmo"),
    ("FUEL_TEMPERATURE", "current_tfu", "base_tfu"),
)
_HISTORY_BRANCHES = (
    ("BORON", "history_bor", "current_bor"),
    ("MODERATOR_TEMPERATURE", "history_tmo", "current_tmo"),
    ("FUEL_TEMPERATURE", "history_tfu", "current_tfu"),
)
_ZERO_SCATTERING_FIELDS = ("sigs_11", "sigs_21", "sigs_22")
_FLOAT_PATTERN = r"[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][-+]?\d+)?"


@dataclass(frozen=True)
class Curve:
    kind: str
    key: float | str | tuple[float, ...] | None
    tables: dict[str, ScalarTable]
    reference_tables: dict[str, ScalarTable] | None = None


@dataclass(frozen=True)
class MaterialLibrary:
    composition_id: int
    composition_name: str
    source_cax: str
    library: CaxLibrary
    fields: list[str]
    curves: list[Curve]
    hm_smeared_density_kg_m3: float | None


def generate_depl_library(
    cax_paths: str | Path | Sequence[str | Path],
    output_path: str | Path,
    *,
    composition_id_start: int = 1,
    name: str | None = None,
    hm_smeared_density_kg_m3: float | None = None,
) -> dict[str, object]:
    """Convert one or more CASMO CAX files into a sparse KOMODO depletion library."""

    if composition_id_start < 1:
        raise ValueError("First composition ID must be a positive integer.")
    if hm_smeared_density_kg_m3 is not None and (
        not math.isfinite(hm_smeared_density_kg_m3) or hm_smeared_density_kg_m3 <= 0.0
    ):
        raise ValueError("Smeared heavy-metal density must be finite and greater than zero.")

    if isinstance(cax_paths, (str, Path)):
        input_paths = [Path(cax_paths)]
    else:
        input_paths = [Path(path) for path in cax_paths]
    if not input_paths:
        raise ValueError("Provide at least one CAX file.")
    if name is not None and len(input_paths) != 1:
        raise ValueError("--name can only be used when converting one CAX file.")

    output_path = Path(output_path)
    output_resolved = output_path.resolve()
    if output_resolved in {path.resolve() for path in input_paths}:
        raise ValueError("Output path must not overwrite one of the input CAX files.")

    material_names = [
        name if name is not None else read_sim_name(cax_path)
        for cax_path in input_paths
    ]
    normalized_names = [material_name.casefold() for material_name in material_names]
    if len(normalized_names) != len(set(normalized_names)):
        raise ValueError("CASMO SIM names must be unique within one depletion library.")

    materials = [
        prepare_material_library(
            cax_path,
            composition_id=composition_id_start + index,
            name=material_names[index],
            hm_smeared_density_kg_m3=hm_smeared_density_kg_m3,
        )
        for index, cax_path in enumerate(input_paths)
    ]
    lines = write_library(materials)
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    return {
        "output_path": str(output_path),
        "materials": [
            {
                "composition_id": material.composition_id,
                "composition_name": material.composition_name,
                "kind": "REFLECTOR" if material.library.is_reflector else "FUEL",
                "base_burnup_points": len(curve_points(material.curves[0], material.fields)),
                "current_curves": sum(curve.kind.startswith("CURRENT_") for curve in material.curves),
                "history_curves": sum(curve.kind.startswith("HISTORY_") for curve in material.curves),
                "rod_curves": sum(curve.kind == "ROD" for curve in material.curves),
                "reflector_interaction_curves": sum(
                    curve.kind == "REFLECTOR_INTERACTION" for curve in material.curves
                ),
                "hm_smeared_density_kg_m3": material.hm_smeared_density_kg_m3,
            }
            for material in materials
        ],
    }


def prepare_material_library(
    cax_path: Path,
    *,
    composition_id: int,
    name: str,
    hm_smeared_density_kg_m3: float | None,
) -> MaterialLibrary:
    if not name.strip():
        raise ValueError(f"Composition name cannot be empty for {cax_path}.")

    states = read_cax_states(cax_path)
    if not states:
        raise ValueError(f"No usable CAX statepoints found in {cax_path}.")
    # SDC TIT records are intentionally excluded from the depletion library.
    states = [state for state in states if "SDC" not in state.title]
    if not states:
        raise ValueError(f"No usable non-SDC CAX statepoints found in {cax_path}.")

    unrodded_states = [state for state in states if state.rod_state == "OUT"]
    if not unrodded_states:
        raise ValueError(f"No unrodded CAX statepoints found in {cax_path}.")

    reflector_interaction_states: list[State] = []
    if unrodded_states[0].is_reflector:
        states, reference_conditions = normalize_reflector_histories(cax_path, states)
        unrodded_states = [state for state in states if state.rod_state == "OUT"]
        require_reflector_reference_output(cax_path, unrodded_states)
        separable_states, reflector_interaction_states = split_reflector_states(
            unrodded_states, reference_conditions=reference_conditions
        )
        library = build_library(separable_states)
    else:
        base_state = unrodded_states[0]
        unclassified = [
            state
            for state in unrodded_states
            if classify_state(
                state,
                base_tfu=base_state.htfu,
                base_tmo=base_state.htmo,
                base_bor=base_state.hbor,
            ) is None
        ]
        if unclassified:
            state = unclassified[0]
            raise ValueError(
                f"Cannot represent non-SDC CAX statepoint at burnup {state.burnup:g} "
                f"with combined branch conditions in {cax_path}: {state.title.strip()}"
            )
        library = build_library(unrodded_states)
    physical_fields = sorted(library.fields)
    fields = sorted([*physical_fields, *_ZERO_SCATTERING_FIELDS])
    base_tables = {field: library.fields[field].base for field in physical_fields}
    base_curve = Curve("BASE", None, base_tables)
    if not curve_points(base_curve, fields):
        raise ValueError(f"The CAX file {cax_path} has no reference burnup curve.")

    curves: list[Curve] = [base_curve]
    curves.extend(make_current_curves(library, physical_fields))
    curves.extend(make_history_curves(library, physical_fields))
    if library.is_reflector:
        if reflector_interaction_states:
            base_history_states = [
                state
                for state in unrodded_states
                if close_conditions(state.htfu, library.base_tfu)
                and close_conditions(state.htmo, library.base_tmo)
                and close_conditions(state.hbor, library.base_bor)
            ]
            if any(
                close_conditions(state.htfu, library.base_tfu)
                and close_conditions(state.htmo, library.base_tmo)
                and close_conditions(state.hbor, library.base_bor)
                for state in reflector_interaction_states
            ):
                # Include separable branch points too; their zero residuals complete
                # the interaction grid for interpolation between stored conditions.
                curves.extend(
                    make_reflector_interaction_curves(
                        base_history_states, library, physical_fields
                    )
                )
            other_interaction_states = [
                state
                for state in reflector_interaction_states
                if not (
                    close_conditions(state.htfu, library.base_tfu)
                    and close_conditions(state.htmo, library.base_tmo)
                    and close_conditions(state.hbor, library.base_bor)
                )
            ]
            curves.extend(
                make_reflector_interaction_curves(
                    other_interaction_states, library, physical_fields
                )
            )
    unsupported_rods = [
        state
        for state in states
        if state.rod_state != "OUT" and not is_direct_rod_state(state, library)
    ]
    if unsupported_rods:
        state = unsupported_rods[0]
        raise ValueError(
            f"Cannot represent rodded CAX statepoint at burnup {state.burnup:g} "
            f"with non-reference branch conditions in {cax_path}: {state.title.strip()}"
        )
    curves.extend(make_rod_curves(states, library, physical_fields))

    if library.is_reflector:
        smeared_density = None
    elif hm_smeared_density_kg_m3 is not None:
        smeared_density = hm_smeared_density_kg_m3
    else:
        cax_lines = cax_path.read_text(encoding="utf-8").splitlines()
        smeared_density = calculate_hm_smeared_density(cax_lines, cax_path)

    return MaterialLibrary(
        composition_id=composition_id,
        composition_name=name,
        source_cax=cax_path.name,
        library=library,
        fields=fields,
        curves=curves,
        hm_smeared_density_kg_m3=smeared_density,
    )


def require_reflector_reference_output(cax_path: Path, reflector_states: list[State]) -> None:
    """Require the first reference statepoint to include reflector-specific XS records."""

    lines = cax_path.read_text(encoding="utf-8").splitlines()
    starts = [index for index, line in enumerate(lines) if line.startswith("TIT")]
    if not starts:
        return

    reference_states = [
        (state.tfu, state.htfu, state.tmo, state.htmo, state.bor, state.hbor)
        for state in reflector_states
    ]
    for position, start in enumerate(starts):
        end = starts[position + 1] if position + 1 < len(starts) else len(lines)
        block = lines[start:end]
        if has_reflector_record(block):
            return
        if len(block) < 3:
            continue

        title = block[0].strip()
        if block[1].strip().startswith("*"):
            title = f"{title} {block[1].strip('* ')}".strip()
        if "SDC" in title:
            continue

        condition = parse_condition_line(block[2])
        if condition is None:
            continue
        _, tfu, htfu, tmo, htmo, bor, hbor = condition
        reference_condition = (tfu, htfu, tmo, htmo, bor, hbor)
        has_matching_reflector_state = any(
            all(
                close_conditions(value, available)
                for value, available in zip(reference_condition, known)
            )
            for known in reference_states
        )
        if has_matching_reflector_state:
            return

        raise ValueError(
            f"Reflector reference statepoint in {cax_path} has no REF/RFH cross-section output "
            f"for TFU/HTFU/TMO/HTMO/BOR/HBOR="
            f"{tfu:g}/{htfu:g}/{tmo:g}/{htmo:g}/{bor:g}/{hbor:g}. "
            "The reference condition must have reflector-specific REF/RFH POL and ADF data; "
            "the ordinary assembly POL/ADF records are not a substitute. "
            "Regenerate the CAX so this reference TIT block contains a REF record with "
            "*POL RFH and *ADF RFH data."
        )


def calculate_hm_smeared_density(cax_lines: list[str], cax_path: Path) -> float:
    """Compute assembly-average HM density from CASMO's HYD line mass and PWR width."""

    hyd_rows = parse_record_rows(cax_lines, "HYD")
    if not hyd_rows or not hyd_rows[0]:
        raise ValueError(f"No usable HYD line-mass value found in {cax_path}.")
    hm_line_mass_g_cm = hyd_rows[0][0]
    if not math.isfinite(hm_line_mass_g_cm) or hm_line_mass_g_cm <= 0.0:
        raise ValueError(f"Invalid HYD line-mass value in {cax_path}.")

    assembly_width_cm = None
    for line in cax_lines:
        parts = line.split()
        if len(parts) < 4 or parts[-1] != "PWR":
            continue
        try:
            assembly_width_cm = float(parts[3])
        except ValueError:
            continue
        break

    if assembly_width_cm is None or not math.isfinite(assembly_width_cm) or assembly_width_cm <= 0.0:
        raise ValueError(
            f"Could not infer square PWR assembly width from {cax_path}; "
            "supply --hm-smeared-density-kg-m3 explicitly."
        )

    # HYD is in g/cm and the PWR cross-sectional assembly area is width squared in cm2.
    # Their ratio is g/cm3; multiply by 1000 to write kg/m3.
    return hm_line_mass_g_cm / (assembly_width_cm * assembly_width_cm) * 1000.0


def read_sim_name(cax_path: Path) -> str:
    """Read the assembly name from the first CASMO SIM output record."""

    lines = cax_path.read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(lines):
        if not line.startswith("SIM") or index + 1 >= len(lines):
            continue
        fields = lines[index + 1].split()
        if fields:
            sim_name = fields[0].strip("'\"")
            if sim_name:
                return sim_name

    raise ValueError(f"No SIM name found in CASMO CAX file {cax_path}.")


def normalize_reflector_histories(
    cax_path: Path,
    states: list[State],
) -> tuple[list[State], tuple[float, float, float] | None]:
    """Recover reflector COE history conditions from the reference point and TIT title."""

    lines = cax_path.read_text(encoding="utf-8").splitlines()
    starts = [index for index, line in enumerate(lines) if line.startswith("TIT")]
    reference_conditions = None
    for position, start in enumerate(starts):
        end = starts[position + 1] if position + 1 < len(starts) else len(lines)
        block = lines[start:end]
        if len(block) < 3:
            continue

        title = block[0].strip()
        if block[1].strip().startswith("*"):
            title = f"{title} {block[1].strip('* ')}".strip()
        if "SDC" in title:
            continue

        condition = parse_condition_line(block[2])
        if condition is None:
            continue
        _, tfu, _, tmo, _, bor, _ = condition
        reference_conditions = (tfu, tmo, bor)
        break

    if reference_conditions is None:
        return states, None

    normalized_states = []
    for state in states:
        if "COE" not in state.title.upper():
            normalized_states.append(state)
            continue

        htfu = title_condition(state.title, "HTFU", reference_conditions[0])
        htmo = title_condition(state.title, "HTMO", reference_conditions[1])
        hbor = title_condition(state.title, "HBOR", reference_conditions[2])
        normalized_states.append(
            replace(state, htfu=htfu, htmo=htmo, hbor=hbor)
        )

    return normalized_states, reference_conditions


def title_condition(title: str, name: str, default: float) -> float:
    match = re.search(
        rf"\b{re.escape(name)}\s*=?\s*({_FLOAT_PATTERN})",
        title,
        flags=re.IGNORECASE,
    )
    return float(match.group(1)) if match is not None else default


def make_current_curves(library: CaxLibrary, fields: list[str]) -> list[Curve]:
    curves: list[Curve] = []
    for label, table_name, base_name in _CURRENT_BRANCHES:
        tables_by_key = branch_tables_by_key(library, fields, table_name)
        base_key = getattr(library, base_name)
        for key, tables in sorted(tables_by_key.items()):
            if key == base_key:
                continue
            curves.append(Curve(f"CURRENT_{label}", key, tables, _base_tables(library, fields)))
    return curves


def make_history_curves(library: CaxLibrary, fields: list[str]) -> list[Curve]:
    curves: list[Curve] = []
    for label, history_name, current_name in _HISTORY_BRANCHES:
        history_maps = branch_tables_by_key(library, fields, history_name)
        current_maps = branch_tables_by_key(library, fields, current_name)
        for key, tables in sorted(history_maps.items()):
            if key not in current_maps:
                raise ValueError(
                    f"Cannot represent {label.lower()} history correction at {key:g}; "
                    f"the matching current-condition branch is absent."
                )
            curves.append(Curve(f"HISTORY_{label}", key, tables, current_maps[key]))
    return curves


def split_reflector_states(
    states: list[State],
    *,
    reference_conditions: tuple[float, float, float] | None = None,
) -> tuple[list[State], list[State]]:
    usable_states = [state for state in states if "SDC" not in state.title]
    if not usable_states:
        raise ValueError("Reflector CAX contains no usable non-SDC statepoints.")

    reference = None
    if reference_conditions is not None:
        reference = next(
            (
                state
                for state in usable_states
                if all(
                    close_conditions(value, expected)
                    for value, expected in zip(
                        (state.tfu, state.tmo, state.bor), reference_conditions
                    )
                )
                and close_conditions(state.tfu, state.htfu)
                and close_conditions(state.tmo, state.htmo)
                and close_conditions(state.bor, state.hbor)
            ),
            None,
        )
    if reference is None:
        reference = next(
            (
                state
                for state in usable_states
                if close_conditions(state.tfu, state.htfu)
                and close_conditions(state.tmo, state.htmo)
                and close_conditions(state.bor, state.hbor)
            ),
            None,
        )
    if reference is None:
        raise ValueError(
            "Reflector CAX has no REF/RFH statepoint matching the reference operating conditions."
        )

    base_conditions = (reference.htfu, reference.htmo, reference.hbor)
    separable_states: list[State] = [reference]
    interaction_states: list[State] = []
    for state in usable_states:
        if state is reference:
            continue
        changed_conditions = 0
        for current, history, base in zip(
            (state.tfu, state.tmo, state.bor),
            (state.htfu, state.htmo, state.hbor),
            base_conditions,
        ):
            if not close_conditions(current, base):
                changed_conditions += 1
            if not close_conditions(history, base):
                changed_conditions += 1
            if changed_conditions > 1:
                break

        if changed_conditions <= 1:
            separable_states.append(state)
        else:
            interaction_states.append(state)

    return separable_states, interaction_states


def make_reflector_interaction_curves(
    states: list[State],
    library: CaxLibrary,
    fields: list[str],
) -> list[Curve]:
    """Store reflector statepoints that the separable base/current curves do not represent."""
    interaction_tables: dict[tuple[float, ...], dict[str, ScalarTable]] = {}

    for state in states:
        key = (state.tfu, state.htfu, state.tmo, state.htmo, state.bor, state.hbor)
        tables = interaction_tables.setdefault(
            key, {field: ScalarTable() for field in fields}
        )
        residuals = {}
        for field in fields:
            matrix = library.fields[field]
            separable_value = matrix.base.at(state.burnup)
            separable_value += branch_delta(
                matrix.current_bor, library.base_bor, state.bor, state.burnup
            )
            separable_value += branch_delta(
                matrix.current_tmo, library.base_tmo, state.tmo, state.burnup
            )
            separable_value += branch_delta(
                matrix.current_tfu, library.base_tfu, state.tfu, state.burnup
            )
            residuals[field] = state.values[field] - separable_value

        if any(state.burnup in tables[field].values for field in fields):
            if all(
                state.burnup in tables[field].values
                and tables[field].values[state.burnup] == residuals[field]
                for field in fields
            ):
                continue
            raise ValueError(
                "Conflicting duplicate reflector statepoint at "
                f"TFU/HTFU/TMO/HTMO/BOR/HBOR={key} and burnup {state.burnup:g}."
            )

        for field in fields:
            tables[field].add(state.burnup, residuals[field])

    return [
        Curve("REFLECTOR_INTERACTION", key, tables)
        for key, tables in sorted(interaction_tables.items())
    ]


def close_conditions(a: float, b: float, *, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol


def branch_tables_by_key(
    library: CaxLibrary,
    fields: list[str],
    table_name: str,
) -> dict[float, dict[str, ScalarTable]]:
    keys = set(getattr(library.fields[fields[0]], table_name))
    for field in fields[1:]:
        field_keys = set(getattr(library.fields[field], table_name))
        if field_keys != keys:
            raise ValueError(f"Branch keys for {table_name} differ for field {field}.")
    return {
        key: {field: getattr(library.fields[field], table_name)[key] for field in fields}
        for key in sorted(keys)
    }


def make_rod_curves(states: list[State], library: CaxLibrary, fields: list[str]) -> list[Curve]:
    rod_tables: dict[str, dict[str, ScalarTable]] = {}
    base_tables = _base_tables(library, fields)

    for state in states:
        if not is_direct_rod_state(state, library):
            continue

        tables = rod_tables.setdefault(state.rod_state, {field: ScalarTable() for field in fields})
        for field in fields:
            table = tables[field]
            if state.burnup in table.values:
                raise ValueError(
                    f"Duplicate direct rod state for {state.rod_state} at burnup {state.burnup:g}."
                )
            reference = lookup_no_extrapolation(
                base_tables[field], state.burnup, f"base {field} for rod state {state.rod_state}"
            )
            table.add(state.burnup, state.values[field] - reference)

    curves = []
    for rod_state, tables in sorted(rod_tables.items()):
        curves.append(Curve("ROD", rod_state, tables))
    return curves


def _base_tables(library: CaxLibrary, fields: list[str]) -> dict[str, ScalarTable]:
    return {field: library.fields[field].base for field in fields}


def curve_points(curve: Curve, fields: list[str]) -> list[tuple[float, list[float]]]:
    missing = [
        field
        for field in fields
        if field not in curve.tables and field not in _ZERO_SCATTERING_FIELDS
    ]
    if missing:
        raise ValueError(f"Curve {curve.kind} is missing fields: {', '.join(missing)}.")

    stored_fields = [field for field in fields if field not in _ZERO_SCATTERING_FIELDS]
    if not stored_fields:
        raise ValueError(f"Curve {curve.kind} does not contain any stored fields.")
    burnup_sets = [set(curve.tables[field].values) for field in stored_fields]
    burnups = sorted(burnup_sets[0])
    if any(points != burnup_sets[0] for points in burnup_sets[1:]):
        raise ValueError(f"Fields in curve {curve.kind} do not share the same burnup points.")

    points: list[tuple[float, list[float]]] = []
    for burnup in burnups:
        values = []
        for field in fields:
            if field in _ZERO_SCATTERING_FIELDS:
                values.append(0.0)
                continue
            value = curve.tables[field].values[burnup]
            if curve.reference_tables is not None:
                reference_table = curve.reference_tables.get(field)
                if reference_table is None:
                    raise ValueError(f"Reference curve for field {field} is missing.")
                value -= lookup_no_extrapolation(
                    reference_table, burnup, f"reference {field} for curve {curve.kind}"
                )
            if not math.isfinite(value):
                raise ValueError(f"Non-finite value for {field} at burnup {burnup:g}.")
            values.append(value)
        points.append((burnup, values))
    return points


def lookup_no_extrapolation(table: ScalarTable, x: float, description: str) -> float:
    values = table.values
    if not values:
        raise ValueError(f"Cannot evaluate empty {description} curve.")
    if x in values:
        return values[x]

    xs = sorted(values)
    if len(xs) < 2 or x < xs[0] or x > xs[-1]:
        raise ValueError(f"Burnup {x:g} is outside the range of {description}.")
    for lower, upper in zip(xs, xs[1:]):
        if lower <= x <= upper:
            y1 = values[lower]
            y2 = values[upper]
            return y1 + (y2 - y1) * (x - lower) / (upper - lower)
    raise ValueError(f"Could not bracket burnup {x:g} in {description}.")


def write_library(materials: list[MaterialLibrary]) -> list[str]:
    lines = [
        "# KOMODO macroscopic depletion library; curve values are ordered by FIELDS.",
        "# CURRENT and HISTORY curves contain additive deltas; ROD curves contain additive deltas.",
        "# HISTORY deltas are history-state values minus the matching current-condition branch.",
        "# REFLECTOR_INTERACTION stores residuals at combined states; keys are TFU HTFU TMO HTMO BOR HBOR.",
        "# Add REFLECTOR_INTERACTION deltas to the matching BASE and CURRENT curves.",
        "KOMODO_DEPL 1",
        "GROUPS 2",
        "BURNUP_UNIT MWd/kgHM",
        "MACRO_XS_UNIT 1/cm",
        "KAPPA_FISSION_UNIT MeV/cm",
        "ADF_UNIT dimensionless",
        "BETA_UNIT dimensionless",
        "LAMBDA_UNIT 1/s",
        "INV_VELOCITY_UNIT s/cm",
        "TEMPERATURE_UNIT K",
        "BORON_UNIT ppm",
        "HM_SMEARED_DENSITY_UNIT kg/m3",
    ]

    for material in materials:
        library = material.library
        fields = material.fields
        lines.extend(
            [
                "",
                f"MATERIAL {material.composition_id} {shlex.quote(material.composition_name)} "
                f"{'REFLECTOR' if library.is_reflector else 'FUEL'}",
                f"SOURCE_CAX {shlex.quote(material.source_cax)}",
                f"REFERENCE_FUEL_TEMPERATURE_K {library.base_tfu:.10g}",
                f"REFERENCE_MODERATOR_TEMPERATURE_K {library.base_tmo:.10g}",
                f"REFERENCE_BORON_PPM {library.base_bor:.10g}",
                (
                    "FISSION_SPECTRUM 0 0"
                    if library.is_reflector
                    else "FISSION_SPECTRUM 1 0"
                ),
                (
                    f"HM_SMEARED_DENSITY_KG_M3 {material.hm_smeared_density_kg_m3:.10g}"
                    if material.hm_smeared_density_kg_m3 is not None
                    else "HM_SMEARED_DENSITY_KG_M3 N/A"
                ),
                f"FIELDS {len(fields)} " + " ".join(fields),
            ]
        )

        for curve in material.curves:
            points = curve_points(curve, fields)
            if curve.kind == "BASE":
                descriptor = "CURVE BASE ABSOLUTE"
            elif curve.kind.startswith("CURRENT_") or curve.kind.startswith("HISTORY_"):
                descriptor = f"CURVE {curve.kind} {float(curve.key):.10g} DELTA"
            elif curve.kind == "REFLECTOR_INTERACTION":
                if not isinstance(curve.key, tuple) or len(curve.key) != 6:
                    raise ValueError("Reflector interaction curve must have six operating-condition keys.")
                key_values = " ".join(f"{value:.10g}" for value in curve.key)
                descriptor = f"CURVE REFLECTOR_INTERACTION TFU_HTFU_TMO_HTMO_BOR_HBOR {key_values} DELTA"
            elif curve.kind == "ROD":
                descriptor = f"CURVE ROD {shlex.quote(str(curve.key))} DELTA"
            else:
                raise ValueError(f"Unsupported curve type {curve.kind}.")

            lines.extend(["", descriptor, f"POINTS {len(points)}"])
            for burnup, values in points:
                formatted_values = " ".join(f"{value:.16E}" for value in values)
                lines.append(f"{burnup:.10g} {formatted_values}")
            lines.append("END_CURVE")

        lines.append("END_MATERIAL")

    lines.append("END_LIBRARY")
    return lines


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Convert CASMO CAX statepoints to a sparse KOMODO depletion library."
    )
    parser.add_argument("cax_paths", nargs="+", help="One or more input CASMO CAX files.")
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        help="Output library path. Defaults to the input .dlib name, or komodo.dlib for multiple files.",
    )
    parser.add_argument(
        "--composition-id-start",
        type=int,
        default=1,
        help="ID for the first composition; IDs increase in input order (default: 1).",
    )
    parser.add_argument(
        "--name",
        dest="composition_name",
        help="Override the CASMO SIM name; only valid with one CAX file.",
    )
    parser.add_argument(
        "--hm-smeared-density-kg-m3",
        "--hm-density-kg-m3",
        type=float,
        dest="hm_smeared_density_kg_m3",
        help="Override assembly-average kgHM/m3; otherwise derive it from CAX HYD and PWR geometry.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    cax_paths = [Path(path) for path in args.cax_paths]
    if args.output_path:
        output_path = Path(args.output_path)
    elif len(cax_paths) == 1:
        output_path = cax_paths[0].with_suffix(".dlib")
    else:
        output_path = Path("komodo.dlib")
    summary = generate_depl_library(
        cax_paths,
        output_path,
        composition_id_start=args.composition_id_start,
        name=args.composition_name,
        hm_smeared_density_kg_m3=args.hm_smeared_density_kg_m3,
    )
    print(summary)


if __name__ == "__main__":
    main()
