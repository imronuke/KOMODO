'''
To convert CASMO CAX files to KOMODO tabular library (tab).
Currently only work for CASMO input with S3C 'HOT'
For reflectors, at least two temperature points must be available.
Try to use for reflectors:

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
from dataclasses import dataclass, replace
from itertools import product
from pathlib import Path
import re


_FLOAT_RE = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?")
_RECORD_RE = re.compile(r"^[A-Z0-9]{3}\b")
_MODERATOR_DENSITY_TABLE = {
    543.15: 0.78106745,
    553.15: 0.76428125,
    563.15: 0.74619716,
    573.15: 0.72650785,
    583.15: 0.70475081,
    593.15: 0.68018488,
    603.15: 0.65150307,
    613.15: 0.61590149,
    617.91: 0.59896404,
}


@dataclass(frozen=True)
class State:
    burnup: float
    tfu: float
    htfu: float
    tmo: float
    htmo: float
    bor: float
    hbor: float
    rod_state: str
    title: str
    is_reflector: bool
    values: dict[str, float]


class ScalarTable:
    def __init__(self) -> None:
        self.values: dict[float, float] = {}

    def add(self, burnup: float, value: float) -> None:
        self.values[burnup] = value

    def at(self, burnup: float) -> float:
        return linear_lookup(self.values, burnup)


class FieldMatrix:
    def __init__(self, base_value: float) -> None:
        self.base_value = base_value
        self.base = ScalarTable()
        self.current_bor: dict[float, ScalarTable] = {}
        self.current_tmo: dict[float, ScalarTable] = {}
        self.current_tfu: dict[float, ScalarTable] = {}
        self.history_bor: dict[float, ScalarTable] = {}
        self.history_tmo: dict[float, ScalarTable] = {}
        self.history_tfu: dict[float, ScalarTable] = {}


@dataclass
class CaxLibrary:
    base_tfu: float
    base_tmo: float
    base_bor: float
    is_reflector: bool
    fields: dict[str, FieldMatrix]
    reflector_grid: dict[str, dict[tuple[float, float, float], ScalarTable]] | None = None
    reflector_grid_incomplete: bool = False

    def burnups(self) -> list[float]:
        any_field = next(iter(self.fields.values()))
        return sorted(any_field.base.values)

    def value_at(self, field: str, burnup: float, *, bor: float, tmo: float, tfu: float) -> float:
        if self.reflector_grid_incomplete:
            raise ValueError(
                "Reflector CAX contains combined operating-condition statepoints, "
                "but they do not form a complete interpolation grid."
            )
        if self.reflector_grid is not None:
            value = interpolate_reflector_grid(
                self.reflector_grid[field], burnup, tfu=tfu, tmo=tmo, bor=bor
            )
            if self.is_reflector:
                matrix = self.fields[field]
                value += history_delta(
                    matrix.history_bor, matrix.current_bor, self.base_bor, bor, burnup
                )
                value += history_delta(
                    matrix.history_tmo, matrix.current_tmo, self.base_tmo, tmo, burnup
                )
                value += history_delta(
                    matrix.history_tfu, matrix.current_tfu, self.base_tfu, tfu, burnup
                )
            return value

        matrix = self.fields[field]
        base = matrix.base.at(burnup)
        return (
            base
            + branch_delta(matrix.current_bor, self.base_bor, bor, burnup)
            + branch_delta(matrix.current_tmo, self.base_tmo, tmo, burnup)
            + branch_delta(matrix.current_tfu, self.base_tfu, tfu, burnup)
            + history_delta(matrix.history_bor, matrix.current_bor, self.base_bor, bor, burnup)
            + history_delta(matrix.history_tmo, matrix.current_tmo, self.base_tmo, tmo, burnup)
            + history_delta(matrix.history_tfu, matrix.current_tfu, self.base_tfu, tfu, burnup)
        )

    def base_value_at(self, field: str, burnup: float) -> float:
        return self.fields[field].base.at(burnup)


@dataclass
class RodAdjustedLibrary:
    base_library: CaxLibrary
    rod_state: str
    rod_deltas: dict[str, ScalarTable]

    @property
    def is_reflector(self) -> bool:
        return self.base_library.is_reflector

    def burnups(self) -> list[float]:
        return self.base_library.burnups()

    def value_at(self, field: str, burnup: float, *, bor: float, tmo: float, tfu: float) -> float:
        return self.base_library.value_at(field, burnup, bor=bor, tmo=tmo, tfu=tfu) + self.rod_deltas[field].at(burnup)

    def base_value_at(self, field: str, burnup: float) -> float:
        return self.base_library.base_value_at(field, burnup) + self.rod_deltas[field].at(burnup)


@dataclass
class CaxBundle:
    unrodded: CaxLibrary
    rodded: RodAdjustedLibrary | None


def generate_tab_library(
    cax_path: str | Path,
    output_path: str | Path,
    *,
    tmo_points: list[float],
    bor_points: list[float],
    tfu_points: list[float],
    burnup_points: list[float] | None = None,
) -> dict[str, object]:
    bundle = read_cax_bundle(cax_path)
    burnups = bundle.unrodded.burnups() if burnup_points is None else burnup_points

    lines: list[str] = []
    write_input_control(lines, bundle.rodded is not None, tmo_points, bor_points, tfu_points)
    for burnup in burnups:
        write_burnup_block(lines, bundle.unrodded, bundle.rodded, burnup, tmo_points, bor_points, tfu_points)

    Path(output_path).write_text("\n".join(lines) + "\n")
    return {
        "output_path": str(output_path),
        "burnups": burnups,
        "fields": sorted(bundle.unrodded.fields),
        "is_reflector": bundle.unrodded.is_reflector,
        "n_tmo": len(tmo_points),
        "n_bor": len(bor_points),
        "n_tfu": len(tfu_points),
        "n_libraries_per_burnup": 2 if bundle.rodded is not None else 1,
        "rodded_state": None if bundle.rodded is None else bundle.rodded.rod_state,
    }


def read_cax_bundle(path: str | Path) -> CaxBundle:
    cax_path = Path(path)
    states = read_cax_states(path)
    if not states:
        raise ValueError(f"No usable CAX statepoints found in {path}.")
    states = [state for state in states if "SDC" not in state.title]
    if not states:
        raise ValueError(f"No usable non-SDC CAX statepoints found in {path}.")

    unrodded_states = [state for state in states if state.rod_state == "OUT"]
    if not unrodded_states:
        raise ValueError(f"No unrodded CAX statepoints found in {path}.")

    reference_conditions = None
    if unrodded_states[0].is_reflector:
        states, reference_conditions = normalize_reflector_histories(cax_path, states)
        unrodded_states = [state for state in states if state.rod_state == "OUT"]
        require_reflector_reference_output(cax_path, unrodded_states)

    states = deduplicate_states(states, cax_path)
    unrodded_states = [state for state in states if state.rod_state == "OUT"]
    if unrodded_states[0].is_reflector:
        reference_state = find_reflector_reference_state(
            unrodded_states, reference_conditions, cax_path
        )
        unrodded_states = [reference_state] + [
            state for state in unrodded_states if state is not reference_state
        ]

    unrodded = build_library(unrodded_states)
    rodded = build_rodded_library(states, unrodded)
    return CaxBundle(unrodded=unrodded, rodded=rodded)


def normalize_reflector_histories(
    cax_path: Path,
    states: list[State],
) -> tuple[list[State], tuple[float, float, float] | None]:
    """Recover reflector history conditions from the reference TIT and COE title."""

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

        normalized_states.append(
            replace(
                state,
                htfu=title_condition(state.title, "HTFU", reference_conditions[0]),
                htmo=title_condition(state.title, "HTMO", reference_conditions[1]),
                hbor=title_condition(state.title, "HBOR", reference_conditions[2]),
            )
        )

    return normalized_states, reference_conditions


def title_condition(title: str, name: str, default: float) -> float:
    match = re.search(
        rf"\b{re.escape(name)}\s*=?\s*({_FLOAT_RE.pattern})",
        title,
        flags=re.IGNORECASE,
    )
    return float(match.group(1)) if match is not None else default


def find_reflector_reference_state(
    states: list[State],
    reference_conditions: tuple[float, float, float] | None,
    cax_path: Path,
) -> State:
    if reference_conditions is None:
        raise ValueError(f"Could not determine reference conditions from {cax_path}.")

    reference = next(
        (
            state
            for state in states
            if all(
                close(value, expected)
                for value, expected in zip(
                    (state.tfu, state.tmo, state.bor), reference_conditions
                )
            )
            and close(state.tfu, state.htfu)
            and close(state.tmo, state.htmo)
            and close(state.bor, state.hbor)
        ),
        None,
    )
    if reference is None:
        tfu, tmo, bor = reference_conditions
        raise ValueError(
            f"Reflector CAX {cax_path} has no REF/RFH statepoint at its reference "
            f"TFU/TMO/BOR={tfu:g}/{tmo:g}/{bor:g}."
        )
    return reference


def require_reflector_reference_output(cax_path: Path, reflector_states: list[State]) -> None:
    """Require REF/RFH cross sections at the CAX reference condition."""

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
        if any(
            all(close(value, available) for value, available in zip(reference_condition, known))
            for known in reference_states
        ):
            return

        raise ValueError(
            f"Reflector reference statepoint in {cax_path} has no REF/RFH cross-section output "
            f"for TFU/HTFU/TMO/HTMO/BOR/HBOR="
            f"{tfu:g}/{htfu:g}/{tmo:g}/{htmo:g}/{bor:g}/{hbor:g}. "
            "The reference condition must have reflector-specific REF/RFH POL and ADF data; "
            "ordinary assembly POL/ADF records are not a substitute. Regenerate the CAX so "
            "the reference TIT condition has REF/RFH output."
        )


def deduplicate_states(states: list[State], cax_path: Path) -> list[State]:
    unique_states: dict[tuple[object, ...], State] = {}
    for state in states:
        key = (
            state.rod_state,
            state.burnup,
            state.tfu,
            state.htfu,
            state.tmo,
            state.htmo,
            state.bor,
            state.hbor,
        )
        previous = unique_states.get(key)
        if previous is None:
            unique_states[key] = state
        elif previous.values != state.values:
            raise ValueError(
                f"Conflicting duplicate CAX statepoint in {cax_path} at "
                f"TFU/HTFU/TMO/HTMO/BOR/HBOR="
                f"{state.tfu:g}/{state.htfu:g}/{state.tmo:g}/{state.htmo:g}/"
                f"{state.bor:g}/{state.hbor:g}, burnup {state.burnup:g}."
            )
    return list(unique_states.values())


def build_library(states: list[State]) -> CaxLibrary:
    base_tfu = states[0].htfu
    base_tmo = states[0].htmo
    base_bor = states[0].hbor
    is_reflector = states[0].is_reflector
    if any(state.is_reflector != is_reflector for state in states):
        raise ValueError("Cannot mix fuel and reflector statepoints in one library.")

    fields = {name: FieldMatrix(base_value=0.0) for name in states[0].values}
    for state in states:
        category = classify_state(state, base_tfu=base_tfu, base_tmo=base_tmo, base_bor=base_bor)
        if category is None:
            history_at_reference = (
                close(state.htfu, base_tfu)
                and close(state.htmo, base_tmo)
                and close(state.hbor, base_bor)
            )
            if is_reflector and "SDC" not in state.title and not history_at_reference:
                raise ValueError(
                    "Cannot represent reflector CAX statepoint with combined current and "
                    f"history conditions in XTAB: {state.title.strip()}"
                )
            continue
        kind, key = category
        for name, value in state.values.items():
            matrix = fields[name]
            if kind == "base":
                matrix.base.add(state.burnup, value)
            elif kind == "current_bor":
                table_for(matrix.current_bor, key).add(state.burnup, value)
            elif kind == "current_tmo":
                table_for(matrix.current_tmo, key).add(state.burnup, value)
            elif kind == "current_tfu":
                table_for(matrix.current_tfu, key).add(state.burnup, value)
            elif kind == "history_bor":
                table_for(matrix.history_bor, key).add(state.burnup, value)
            elif kind == "history_tmo":
                table_for(matrix.history_tmo, key).add(state.burnup, value)
            elif kind == "history_tfu":
                table_for(matrix.history_tfu, key).add(state.burnup, value)

    for matrix in fields.values():
        for history_tables, current_tables, name in (
            (matrix.history_bor, matrix.current_bor, "boron"),
            (matrix.history_tmo, matrix.current_tmo, "moderator temperature"),
            (matrix.history_tfu, matrix.current_tfu, "fuel temperature"),
        ):
            missing_current = set(history_tables) - set(current_tables)
            if missing_current:
                key = min(missing_current)
                raise ValueError(
                    f"Cannot represent {name} history correction at {key:g}; "
                    "the matching current-condition branch is absent."
                )
        matrix.current_bor[base_bor] = matrix.base
        matrix.current_tmo[base_tmo] = matrix.base
        matrix.current_tfu[base_tfu] = matrix.base

    reflector_grid = None
    reflector_grid_incomplete = False
    if is_reflector:
        reflector_grid, reflector_grid_incomplete = build_reflector_grid(states)

    return CaxLibrary(
        base_tfu=base_tfu,
        base_tmo=base_tmo,
        base_bor=base_bor,
        is_reflector=is_reflector,
        fields=fields,
        reflector_grid=reflector_grid,
        reflector_grid_incomplete=reflector_grid_incomplete,
    )


def build_reflector_grid(
    states: list[State],
) -> tuple[dict[str, dict[tuple[float, float, float], ScalarTable]] | None, bool]:
    if not states:
        return None, False

    base_tfu, base_tmo, base_bor = states[0].htfu, states[0].htmo, states[0].hbor
    current_states = [
        state
        for state in states
        if "SDC" not in state.title
        and close(state.htfu, base_tfu)
        and close(state.htmo, base_tmo)
        and close(state.hbor, base_bor)
    ]
    if not current_states:
        return None, False

    coordinates = {
        (state.tfu, state.tmo, state.bor)
        for state in current_states
    }
    axes = [sorted({coordinate[index] for coordinate in coordinates}) for index in range(3)]
    expected_coordinates = set(product(*axes))
    has_combined_states = any(
        sum(
            not close(value, base)
            for value, base in zip(
                (state.tfu, state.tmo, state.bor),
                (base_tfu, base_tmo, base_bor),
            )
        ) > 1
        for state in current_states
    )
    if coordinates != expected_coordinates:
        return None, has_combined_states

    grid = {name: {} for name in current_states[0].values}
    for state in current_states:
        coordinate = (state.tfu, state.tmo, state.bor)
        for name, value in state.values.items():
            table = grid[name].setdefault(coordinate, ScalarTable())
            if state.burnup in table.values and table.values[state.burnup] != value:
                raise ValueError(
                    "Conflicting duplicate reflector statepoint at "
                    f"TFU/TMO/BOR={state.tfu:g}/{state.tmo:g}/{state.bor:g}, "
                    f"burnup {state.burnup:g}, field {name}."
                )
            table.add(state.burnup, value)

    return grid, False


def build_rodded_library(states: list[State], unrodded: CaxLibrary) -> RodAdjustedLibrary | None:
    rod_states = [state for state in states if is_direct_rod_state(state, unrodded)]
    if not rod_states:
        return None

    rod_labels = sorted({state.rod_state for state in rod_states})
    if len(rod_labels) != 1:
        raise ValueError(f"Expected one direct rod branch, found {rod_labels}.")

    rod_deltas = {name: ScalarTable() for name in unrodded.fields}
    for state in rod_states:
        for name, value in state.values.items():
            delta = value - unrodded.base_value_at(name, state.burnup)
            rod_deltas[name].add(state.burnup, delta)

    return RodAdjustedLibrary(base_library=unrodded, rod_state=rod_labels[0], rod_deltas=rod_deltas)


def read_cax_states(path: str | Path) -> list[State]:
    lines = Path(path).read_text().splitlines()
    starts = [i for i, line in enumerate(lines) if line.startswith("TIT")]
    rod_contexts = rod_context_by_start(lines)
    has_reflector_states = any(line.startswith("REF") for line in lines)
    states: list[State] = []

    for pos, start in enumerate(starts):
        end = starts[pos + 1] if pos + 1 < len(starts) else len(lines)
        block = lines[start:end]
        state = parse_state_block(block, rod_state=rod_contexts.get(start, "OUT"), prefer_reflector=has_reflector_states)
        if state is not None:
            states.append(state)

    return states


def rod_context_by_start(lines: list[str]) -> dict[int, str]:
    contexts: dict[int, str] = {}
    current_rod_state = "OUT"

    for index, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("*I TTL") and " ROD" not in stripped:
            current_rod_state = "OUT"
        elif stripped.startswith("*I ROD "):
            match = re.search(r"'([^']+)'", stripped)
            if match is not None:
                label = match.group(1)
                current_rod_state = "OUT" if label == "OUT" else label
        if line.startswith("TIT"):
            contexts[index] = current_rod_state

    return contexts


def parse_state_block(block: list[str], *, rod_state: str, prefer_reflector: bool) -> State | None:
    if len(block) < 4:
        return None

    title = block[0].strip()
    if block[1].strip().startswith("*"):
        title = f"{title} {block[1].strip('* ')}".strip()

    condition = parse_condition_line(block[2])
    if condition is None:
        return None
    burnup, tfu, htfu, tmo, htmo, bor, hbor = condition

    is_reflector = has_reflector_record(block)
    if prefer_reflector and not is_reflector:
        return None

    if is_reflector:
        pol = parse_reflector_record_rows(block, "POL", "RFH")
        adf = parse_reflector_record_rows(block, "ADF", "RFH")
        if not pol or not adf:
            raise ValueError("Reflector CAX state does not contain REF/RFH POL and ADF records.")
        values = parse_reflector_pol_values(pol)
        values.update(parse_reflector_adf_values(adf))
        values.update(zero_delayed_values())
    else:
        pol = parse_record_rows(block, "POL")
        adf = parse_record_rows(block, "ADF")
        delayed = parse_record_rows(block, "DEL")
        if not pol or not adf or not delayed:
            return None

        values = parse_pol_values(pol)
        values.update(parse_adf_values(adf))
        values.update(parse_del_values(delayed))

    return State(
        burnup=burnup,
        tfu=tfu,
        htfu=htfu,
        tmo=tmo,
        htmo=htmo,
        bor=bor,
        hbor=hbor,
        rod_state=rod_state,
        title=title,
        is_reflector=is_reflector,
        values=values,
    )


def parse_condition_line(line: str) -> tuple[float, float, float, float, float, float, float] | None:
    try:
        burnup = float(line.split()[0])
    except ValueError:
        return None
    pairs = re.findall(r"([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?)\s*/\s*([-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[Ee][-+]?\d+)?)", line)
    if len(pairs) < 4:
        return None
    tfu, htfu = float(pairs[1][0]), float(pairs[1][1])
    tmo, htmo = float(pairs[2][0]), float(pairs[2][1])
    bor, hbor = float(pairs[3][0]), float(pairs[3][1])
    return burnup, tfu, htfu, tmo, htmo, bor, hbor


def parse_pair(token: str) -> tuple[float, float]:
    left, right = token.split("/", 1)
    return float(left), float(right)


def parse_record_rows(block: list[str], name: str) -> list[list[float]]:
    for index, line in enumerate(block):
        if not line.startswith(name):
            continue
        rows: list[list[float]] = []
        for data_line in block[index + 1 :]:
            if _RECORD_RE.match(data_line):
                break
            numbers = floats_in(data_line)
            if numbers:
                rows.append(numbers)
        return rows
    return []


def has_reflector_record(block: list[str]) -> bool:
    return any(line.startswith("REF") for line in block)


def parse_reflector_record_rows(block: list[str], name: str, region: str) -> list[list[float]]:
    in_ref = False
    for index, line in enumerate(block):
        if line.startswith("REF"):
            in_ref = True
            continue
        if not in_ref:
            continue
        if _RECORD_RE.match(line):
            break

        parts = line.strip().split()
        if len(parts) < 2 or parts[0] != f"*{name}" or parts[1] != region:
            continue

        rows: list[list[float]] = []
        for data_line in block[index + 1 :]:
            if data_line.startswith("*") or _RECORD_RE.match(data_line):
                break
            numbers = floats_in(data_line)
            if numbers:
                rows.append(numbers)
        return rows
    return []


def parse_pol_values(rows: list[list[float]]) -> dict[str, float]:
    flat = [value for row in rows for value in row]
    if len(flat) < 34:
        raise ValueError("POL record does not contain the expected two-group constants.")

    d1, d2, casmo_remov1, siga1, siga2, nusigf1, nusigf2 = flat[:7]
    nu1, nu2, kappa1, kappa2 = flat[-4:]
    if nu1 == 0.0 or nu2 == 0.0:
        raise ValueError("POL record contains zero group-wise nu; cannot compute kappa-fission.")

    # CASMO POL labels this value REMOV1; in the two-group balance it is the
    # group-1 to group-2 removal/downscatter term.
    sigs12 = casmo_remov1
    return {
        "transport_1": 1.0 / (3.0 * d1),
        "transport_2": 1.0 / (3.0 * d2),
        "siga_1": siga1,
        "siga_2": siga2,
        "nusigf_1": nusigf1,
        "nusigf_2": nusigf2,
        "kappa_fission_1": kappa1 * nusigf1 / nu1,
        "kappa_fission_2": kappa2 * nusigf2 / nu2,
        "sigs_12": sigs12,
    }


def parse_reflector_pol_values(rows: list[list[float]]) -> dict[str, float]:
    flat = [value for row in rows for value in row]
    if len(flat) < 11:
        raise ValueError("REF/RFH POL record does not contain the expected two-group constants.")

    d1, d2, casmo_remov1, siga1, siga2, nusigf1, nusigf2 = flat[:7]
    tail = flat[7:11]
    return {
        "transport_1": 1.0 / (3.0 * d1),
        "transport_2": 1.0 / (3.0 * d2),
        "siga_1": siga1,
        "siga_2": siga2,
        "nusigf_1": nusigf1,
        "nusigf_2": nusigf2,
        "kappa_fission_1": 0.0,
        "kappa_fission_2": 0.0,
        "sigs_12": casmo_remov1,
        "inv_velocity_1": tail[1],
        "inv_velocity_2": tail[3],
    }


def parse_adf_values(rows: list[list[float]]) -> dict[str, float]:
    adf_rows = [row for row in rows if len(row) >= 5]
    if len(adf_rows) < 2:
        raise ValueError("ADF record does not contain ADF1 and ADF2 rows.")
    return {
        "adf_1": adf_rows[0][4],
        "adf_2": adf_rows[1][4],
    }


def parse_reflector_adf_values(rows: list[list[float]]) -> dict[str, float]:
    adf_rows = [row for row in rows if len(row) >= 1]
    if len(adf_rows) < 2:
        raise ValueError("REF/RFH ADF record does not contain ADF1 and ADF2 rows.")
    return {
        "adf_1": adf_rows[0][0],
        "adf_2": adf_rows[1][0],
    }


def parse_del_values(rows: list[list[float]]) -> dict[str, float]:
    if len(rows) < 5 or len(rows[1]) < 6 or len(rows[2]) < 6 or len(rows[4]) < 2:
        raise ValueError("DEL record does not contain beta, lambda, and inverse velocity rows.")

    values: dict[str, float] = {
        "inv_velocity_1": rows[4][0],
        "inv_velocity_2": rows[4][1],
    }
    for index, value in enumerate(rows[2][:6], start=1):
        values[f"lambda_{index}"] = value
    for index, value in enumerate(rows[1][:6], start=1):
        values[f"beta_{index}"] = value
    return values


def zero_delayed_values() -> dict[str, float]:
    values: dict[str, float] = {}
    for index in range(1, 7):
        values[f"lambda_{index}"] = 0.0
        values[f"beta_{index}"] = 0.0
    return values


def classify_state(
    state: State,
    *,
    base_tfu: float,
    base_tmo: float,
    base_bor: float,
) -> tuple[str, float | None] | None:
    if "SDC" in state.title:
        return None

    tfu_base = close(state.tfu, base_tfu) and close(state.htfu, base_tfu)
    tmo_base = close(state.tmo, base_tmo) and close(state.htmo, base_tmo)
    bor_base = close(state.bor, base_bor) and close(state.hbor, base_bor)

    if tfu_base and tmo_base and bor_base:
        return "base", None
    if state.is_reflector:
        if tfu_base and tmo_base and close(state.hbor, base_bor) and not close(state.bor, base_bor):
            return "current_bor", state.bor
        if tfu_base and bor_base and close(state.htmo, base_tmo) and not close(state.tmo, base_tmo):
            return "current_tmo", state.tmo
        if tmo_base and bor_base and close(state.htfu, base_tfu) and not close(state.tfu, base_tfu):
            return "current_tfu", state.tfu
        if tfu_base and tmo_base and close(state.bor, state.hbor) and not close(state.hbor, base_bor):
            return "history_bor", state.hbor
        if tfu_base and bor_base and close(state.tmo, state.htmo) and not close(state.htmo, base_tmo):
            return "history_tmo", state.htmo
        if tmo_base and bor_base and close(state.tfu, state.htfu) and not close(state.htfu, base_tfu):
            return "history_tfu", state.htfu
        return None
    if tfu_base and tmo_base and close(state.hbor, base_bor) and not close(state.bor, base_bor):
        return "current_bor", state.bor
    if tfu_base and bor_base and close(state.htmo, base_tmo) and not close(state.tmo, base_tmo):
        return "current_tmo", state.tmo
    if tmo_base and bor_base and close(state.htfu, base_tfu) and not close(state.tfu, base_tfu):
        return "current_tfu", state.tfu
    if tfu_base and tmo_base and close(state.bor, state.hbor) and not close(state.hbor, base_bor):
        return "history_bor", state.hbor
    if tfu_base and bor_base and close(state.tmo, state.htmo) and not close(state.htmo, base_tmo):
        return "history_tmo", state.htmo
    if tmo_base and bor_base and close(state.tfu, state.htfu) and not close(state.htfu, base_tfu):
        return "history_tfu", state.htfu
    return None


def is_direct_rod_state(state: State, library: CaxLibrary) -> bool:
    if state.rod_state == "OUT":
        return False

    tfu_base = close(state.tfu, library.base_tfu) and close(state.htfu, library.base_tfu)
    tmo_base = close(state.tmo, library.base_tmo) and close(state.htmo, library.base_tmo)
    bor_base = close(state.bor, library.base_bor) and close(state.hbor, library.base_bor)
    return tfu_base and tmo_base and bor_base and "ROD" in state.title


def write_burnup_block(
    lines: list[str],
    library: CaxLibrary,
    rodded_library: RodAdjustedLibrary | None,
    burnup: float,
    tmo_points: list[float],
    bor_points: list[float],
    tfu_points: list[float],
) -> None:
    lines.extend(
        [
            "",
            " * --------------------------------------------------------",
            f" * BURNUP {burnup:6.2f}",
            " * --------------------------------------------------------",
            " *",
        ]
    )
    write_library_tables(lines, library, burnup, tmo_points, bor_points, tfu_points)
    if rodded_library is not None:
        write_library_tables(lines, rodded_library, burnup, tmo_points, bor_points, tfu_points)

    vector_library = library if rodded_library is None else rodded_library
    fission_spectrum = [0.0, 0.0] if library.is_reflector else [1.0, 0.0]
    write_constant_vector(lines, "Fission Spectrum", fission_spectrum, groups=2)
    write_constant_vector(
        lines,
        "Inverse velocity",
        [vector_library.base_value_at("inv_velocity_1", burnup), vector_library.base_value_at("inv_velocity_2", burnup)],
        groups=2,
    )
    write_constant_vector(
        lines,
        "Delay Neutron Decay Constant (Lambda)",
        [vector_library.base_value_at(f"lambda_{i}", burnup) for i in range(1, 7)],
        groups=6,
    )
    write_constant_vector(
        lines,
        "Delay Neutron Fraction (Beta)",
        [vector_library.base_value_at(f"beta_{i}", burnup) for i in range(1, 7)],
        groups=6,
    )


def write_input_control(
    lines: list[str],
    has_control_rods: bool,
    tmo_points: list[float],
    bor_points: list[float],
    tfu_points: list[float],
) -> None:
    mod_densities = [moderator_density_from_temperature(tmo) for tmo in tmo_points]
    cr_flag = 1 if has_control_rods else 0
    lines.extend(
        [
            "  * INPUT CONTROL",
            f"   1   {cr_flag} * ADF and CR",
            " *     Mod Dens      Boron ppm      Fuel Temp       Mod Temp",
            f"{len(tmo_points):12d}{len(bor_points):13d}{len(tfu_points):13d}{1:13d}",
            format_values_fixed(mod_densities),
            format_values_fixed(bor_points),
            format_values_fixed(tfu_points),
        ]
    )


def write_library_tables(
    lines: list[str],
    library,
    burnup: float,
    tmo_points: list[float],
    bor_points: list[float],
    tfu_points: list[float],
) -> None:
    write_group_table(lines, "Transport XSEC Table", [(1, "transport_1"), (2, "transport_2")], library, burnup, tmo_points, bor_points, tfu_points)
    write_group_table(lines, "Absorption XSEC Table", [(1, "siga_1"), (2, "siga_2")], library, burnup, tmo_points, bor_points, tfu_points)
    write_group_table(lines, "NU-Fission XSEC Table", [(1, "nusigf_1"), (2, "nusigf_2")], library, burnup, tmo_points, bor_points, tfu_points)
    write_group_table(lines, "Kappa-Fission XSEC Table", [(1, "kappa_fission_1"), (2, "kappa_fission_2")], library, burnup, tmo_points, bor_points, tfu_points)
    write_scattering_table(lines, library, burnup, tmo_points, bor_points, tfu_points)
    write_group_table(lines, "ADF Table", [(1, "adf_1"), (2, "adf_2")], library, burnup, tmo_points, bor_points, tfu_points)


def write_group_table(
    lines: list[str],
    title: str,
    groups: list[tuple[int, str]],
    library: CaxLibrary,
    burnup: float,
    tmo_points: list[float],
    bor_points: list[float],
    tfu_points: list[float],
) -> None:
    lines.extend([" *", f" * {title}", " *"])
    for group, field in groups:
        lines.append(f"  * GROUP {group:6d}")
        write_grid(lines, lambda tfu, bor, tmo, field=field: library.value_at(field, burnup, bor=bor, tmo=tmo, tfu=tfu), tmo_points, bor_points, tfu_points)


def write_zero_group_table(
    lines: list[str],
    title: str,
    groups: list[int],
    tmo_points: list[float],
    bor_points: list[float],
    tfu_points: list[float],
) -> None:
    lines.extend([" *", f" * {title}", " *"])
    for group in groups:
        lines.append(f"  * GROUP {group:6d}")
        write_grid(lines, lambda _tfu, _bor, _tmo: 0.0, tmo_points, bor_points, tfu_points)


def write_scattering_table(
    lines: list[str],
    library: CaxLibrary,
    burnup: float,
    tmo_points: list[float],
    bor_points: list[float],
    tfu_points: list[float],
) -> None:
    lines.extend([" *", " * Scattering XSEC Table", " *"])
    for from_group, to_group, field in [(1, 1, None), (1, 2, "sigs_12"), (2, 1, None), (2, 2, None)]:
        lines.append(f"  * GROUP {from_group:6d} -> {to_group:4d}")
        if field is None:
            write_grid(lines, lambda _tfu, _bor, _tmo: 0.0, tmo_points, bor_points, tfu_points)
        else:
            write_grid(lines, lambda tfu, bor, tmo, field=field: library.value_at(field, burnup, bor=bor, tmo=tmo, tfu=tfu), tmo_points, bor_points, tfu_points)


def write_grid(lines: list[str], value_fn, tmo_points: list[float], bor_points: list[float], tfu_points: list[float]) -> None:
    for tfu in tfu_points:
        for bor in bor_points:
            values = [value_fn(tfu, bor, tmo) for tmo in tmo_points]
            lines.append(format_values(values))


def write_constant_vector(lines: list[str], title: str, values: list[float], *, groups: int) -> None:
    lines.extend([" *", f" * {title}", " *"])
    header = "".join(f"{group:16d}" for group in range(1, groups + 1))
    lines.append(f"  * GROUP {header}")
    lines.append(format_values(values))


def format_values(values: list[float]) -> str:
    return "".join(f"{value:16.5E}" for value in values)


def format_values_fixed(values: list[float]) -> str:
    return "".join(f"{value:16.5f}" for value in values)


def floats_in(line: str) -> list[float]:
    return [float(value) for value in _FLOAT_RE.findall(line)]


def table_for(tables: dict[float, ScalarTable], key: float | None) -> ScalarTable:
    if key is None:
        raise ValueError("Cannot create branch table without a branch key.")
    if key not in tables:
        tables[key] = ScalarTable()
    return tables[key]


def close(a: float, b: float, *, tol: float = 1e-6) -> bool:
    return abs(a - b) <= tol


def linear_lookup(values: dict[float, float], x: float) -> float:
    if not values:
        raise ValueError("Cannot interpolate an empty table.")
    if x in values:
        return values[x]
    xs = sorted(values)
    if len(xs) == 1:
        return values[xs[0]]
    if x < xs[0]:
        x1, x2 = xs[0], xs[1]
    elif x > xs[-1]:
        x1, x2 = xs[-2], xs[-1]
    else:
        for lower, upper in zip(xs, xs[1:]):
            if lower <= x <= upper:
                x1, x2 = lower, upper
                break
        else:
            raise ValueError(f"Could not bracket {x}.")
    y1 = values[x1]
    y2 = values[x2]
    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)


def branch_delta(tables: dict[float, ScalarTable], base_value: float, desired_value: float, burnup: float) -> float:
    branch_value = interpolate_tables(tables, desired_value, burnup)
    return branch_value - tables[base_value].at(burnup)


def history_delta(
    history_tables: dict[float, ScalarTable],
    branch_tables: dict[float, ScalarTable],
    base_value: float,
    desired_value: float,
    burnup: float,
) -> float:
    corrections = {base_value: 0.0}
    for history_value, history_table in history_tables.items():
        if history_value in branch_tables:
            corrections[history_value] = history_table.at(burnup) - branch_tables[history_value].at(burnup)
    return linear_lookup(corrections, desired_value)


def interpolate_tables(tables: dict[float, ScalarTable], desired_value: float, burnup: float) -> float:
    values = {key: table.at(burnup) for key, table in tables.items()}
    return linear_lookup(values, desired_value)


def interpolate_reflector_grid(
    tables: dict[tuple[float, float, float], ScalarTable],
    burnup: float,
    *,
    tfu: float,
    tmo: float,
    bor: float,
) -> float:
    coordinates = list(tables)
    axes = [sorted({coordinate[index] for coordinate in coordinates}) for index in range(3)]
    target = (tfu, tmo, bor)
    weights_by_axis = [linear_weights(axis, value) for axis, value in zip(axes, target)]

    result = 0.0
    for corner in product(*weights_by_axis):
        coordinate = tuple(point for point, _ in corner)
        weight = 1.0
        for _, axis_weight in corner:
            weight *= axis_weight
        result += weight * tables[coordinate].at(burnup)
    return result


def linear_weights(values: list[float], x: float) -> list[tuple[float, float]]:
    if not values:
        raise ValueError("Cannot interpolate an empty grid axis.")
    if len(values) == 1:
        return [(values[0], 1.0)]
    if x in values:
        return [(x, 1.0)]

    if x < values[0]:
        lower, upper = values[:2]
    elif x > values[-1]:
        lower, upper = values[-2:]
    else:
        for lower, upper in zip(values, values[1:]):
            if lower < x < upper:
                break
        else:
            raise ValueError(f"Could not bracket {x} on reflector grid axis.")

    upper_weight = (x - lower) / (upper - lower)
    return [(lower, 1.0 - upper_weight), (upper, upper_weight)]


def moderator_density_from_temperature(temperature_k: float) -> float:
    return linear_lookup(_MODERATOR_DENSITY_TABLE, temperature_k)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate a .tab MGXS file from a CASMO .cax file.")
    parser.add_argument("cax_path", help="Path to the CASMO CAX file.")
    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        help="Output TAB file path. Defaults to the input filename with a .tab extension.",
    )
    return parser.parse_args()


def default_output_path(cax_path: str | Path) -> Path:
    return Path(cax_path).with_suffix(".tab")


def main() -> None:
    args = parse_args()
    cax_path = Path(args.cax_path)
    output_path = Path(args.output_path) if args.output_path else default_output_path(cax_path)

    summary = generate_tab_library(
        cax_path,
        output_path,
        tmo_points=[600.0, 580.0, 560.0],
        bor_points=[0.0, 1000.0, 2400.0],
        tfu_points=[600.0, 900.0, 1200.0],
        burnup_points=[0.0, 5.0, 50.0],
    )
    print(summary)


if __name__ == "__main__":
    main()
