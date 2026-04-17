#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Usage:
  baseline_stage0.sh run
  baseline_stage0.sh check

Environment overrides:
  KOMODO_BIN                     Path to komodo executable (default: build/komodo)
  BASELINE_STAGE0_ARTIFACT_DIR   Output directory (default: /tmp/komodo_stage0_XXXXXXXX)
  BASELINE_STAGE0_CASES          Comma-separated subset (IAEA3Ds,A1t_coarse,A1,A2,adjoint,BIBLIS,PNM), default all

  BASELINE_STAGE0_EXPECT_IAEA3DS_KEFF
  BASELINE_STAGE0_EXPECT_A1_BORON
  BASELINE_STAGE0_EXPECT_A1_KEFF
  BASELINE_STAGE0_EXPECT_A1T_COARSE_PEAK_STEP
  BASELINE_STAGE0_EXPECT_A1T_COARSE_PEAK_REL_POWER
  BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_STEP
  BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_TIME
  BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_REACT
  BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_REL_POWER
  BASELINE_STAGE0_EXPECT_A1T_COARSE_AVG_FUEL_TEMP
  BASELINE_STAGE0_EXPECT_A1T_COARSE_MAX_FUEL_TEMP
  BASELINE_STAGE0_EXPECT_A1T_COARSE_AVG_MOD_TEMP
  BASELINE_STAGE0_EXPECT_A1T_COARSE_MAX_MOD_TEMP
  BASELINE_STAGE0_EXPECT_A1T_COARSE_OUTLET_MOD_TEMP
  BASELINE_STAGE0_EXPECT_A1T_COARSE_AVG_MOD_DENS
  BASELINE_STAGE0_EXPECT_A1T_COARSE_OUTLET_MOD_DENS
  BASELINE_STAGE0_EXPECT_A2_BORON
  BASELINE_STAGE0_EXPECT_A2_KEFF
  BASELINE_STAGE0_EXPECT_ADJOINT_KEFF
  BASELINE_STAGE0_EXPECT_BIBLIS_KEFF
  BASELINE_STAGE0_EXPECT_PNM_KEFF

  BASELINE_STAGE0_TOL_KEFF
  BASELINE_STAGE0_TOL_BORON
  BASELINE_STAGE0_TOL_REACTIVITY
  BASELINE_STAGE0_TOL_REL_POWER
EOF
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT_DIR="${SCRIPT_DIR}"
KOMODO_BIN="${KOMODO_BIN:-${ROOT_DIR}/build/komodo}"

if [[ ! -x "${KOMODO_BIN}" ]]; then
  echo "ERROR: KOMODO executable not found or not executable: ${KOMODO_BIN}" >&2
  exit 2
fi

if [[ -n "${BASELINE_STAGE0_ARTIFACT_DIR:-}" ]]; then
  ARTIFACT_DIR="${BASELINE_STAGE0_ARTIFACT_DIR}"
else
  ARTIFACT_DIR="$(mktemp -d /tmp/komodo_stage0_XXXXXXXX)"
fi

RUN_DIR="${ARTIFACT_DIR}/run"
METRICS_FILE="${ARTIFACT_DIR}/metrics.tsv"
CASES="${BASELINE_STAGE0_CASES:-IAEA3Ds,A1t_coarse,A1,A2,adjoint,BIBLIS,PNM}"

expect_iaea3ds_keff="${BASELINE_STAGE0_EXPECT_IAEA3DS_KEFF:-1.029084}"
expect_a1_boron="${BASELINE_STAGE0_EXPECT_A1_BORON:-560.52}"
expect_a1_keff="${BASELINE_STAGE0_EXPECT_A1_KEFF:-1.00000}"
expect_a1t_coarse_peak_step="${BASELINE_STAGE0_EXPECT_A1T_COARSE_PEAK_STEP:-134}"
expect_a1t_coarse_peak_rel_power="${BASELINE_STAGE0_EXPECT_A1T_COARSE_PEAK_REL_POWER:-7.9681E-01}"
expect_a1t_coarse_final_step="${BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_STEP:-280}"
expect_a1t_coarse_final_time="${BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_TIME:-5.000}"
expect_a1t_coarse_final_react="${BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_REACT:-0.3582}"
expect_a1t_coarse_final_rel_power="${BASELINE_STAGE0_EXPECT_A1T_COARSE_FINAL_REL_POWER:-1.9556E-01}"
expect_a1t_coarse_avg_fuel_temp="${BASELINE_STAGE0_EXPECT_A1T_COARSE_AVG_FUEL_TEMP:-597.0}"
expect_a1t_coarse_max_fuel_temp="${BASELINE_STAGE0_EXPECT_A1T_COARSE_MAX_FUEL_TEMP:-939.0}"
expect_a1t_coarse_avg_mod_temp="${BASELINE_STAGE0_EXPECT_A1T_COARSE_AVG_MOD_TEMP:-563.0}"
expect_a1t_coarse_max_mod_temp="${BASELINE_STAGE0_EXPECT_A1T_COARSE_MAX_MOD_TEMP:-581.3}"
expect_a1t_coarse_outlet_mod_temp="${BASELINE_STAGE0_EXPECT_A1T_COARSE_OUTLET_MOD_TEMP:-566.4}"
expect_a1t_coarse_avg_mod_dens="${BASELINE_STAGE0_EXPECT_A1T_COARSE_AVG_MOD_DENS:-746.2}"
expect_a1t_coarse_outlet_mod_dens="${BASELINE_STAGE0_EXPECT_A1T_COARSE_OUTLET_MOD_DENS:-739.8}"
expect_a2_boron="${BASELINE_STAGE0_EXPECT_A2_BORON:-1156.13}"
expect_a2_keff="${BASELINE_STAGE0_EXPECT_A2_KEFF:-1.00000}"
expect_adjoint_keff="${BASELINE_STAGE0_EXPECT_ADJOINT_KEFF:-1.029090}"
expect_biblis_keff="${BASELINE_STAGE0_EXPECT_BIBLIS_KEFF:-1.025119}"
expect_pnm_keff="${BASELINE_STAGE0_EXPECT_PNM_KEFF:-1.029604}"

tol_keff="${BASELINE_STAGE0_TOL_KEFF:-5e-6}"
tol_boron="${BASELINE_STAGE0_TOL_BORON:-0.10}"
tol_reactivity="${BASELINE_STAGE0_TOL_REACTIVITY:-5e-4}"
tol_rel_power="${BASELINE_STAGE0_TOL_REL_POWER:-1e-2}"

run_stage0() {
  mkdir -p "${RUN_DIR}"

  if [[ "${CASES}" == *"IAEA3Ds"* ]]; then
    cp "${ROOT_DIR}/smpl/static/IAEA3Ds" "${RUN_DIR}/IAEA3Ds"
    "${KOMODO_BIN}" "${RUN_DIR}/IAEA3Ds" > "${RUN_DIR}/IAEA3Ds.stdout" 2>&1
  fi

  if [[ "${CASES}" == *"A1t_coarse"* ]]; then
    cp "${ROOT_DIR}/smpl/transient/NEACRP/A1t_coarse" "${RUN_DIR}/A1t_coarse"
    "${KOMODO_BIN}" "${RUN_DIR}/A1t_coarse" > "${RUN_DIR}/A1t_coarse.stdout" 2>&1
  fi

  if [[ "${CASES}" == *"A1"* ]]; then
    sed 's#/home/imronuke/github/tsuiyong/KOMODO#/home/imronuke/git/KOMODO#g' \
      "${ROOT_DIR}/smpl/static/NEACRP/A1" > "${RUN_DIR}/A1"
    "${KOMODO_BIN}" "${RUN_DIR}/A1" > "${RUN_DIR}/A1.stdout" 2>&1
  fi

  if [[ "${CASES}" == *"A2"* ]]; then
    cp "${ROOT_DIR}/smpl/static/NEACRP/A2" "${RUN_DIR}/A2"
    "${KOMODO_BIN}" "${RUN_DIR}/A2" > "${RUN_DIR}/A2.stdout" 2>&1
  fi

  if [[ "${CASES}" == *"adjoint"* ]]; then
    cp "${ROOT_DIR}/smpl/static/adjoint" "${RUN_DIR}/adjoint"
    "${KOMODO_BIN}" "${RUN_DIR}/adjoint" > "${RUN_DIR}/adjoint.stdout" 2>&1
  fi

  if [[ "${CASES}" == *"BIBLIS"* ]]; then
    cp "${ROOT_DIR}/smpl/static/BIBLIS" "${RUN_DIR}/BIBLIS"
    "${KOMODO_BIN}" "${RUN_DIR}/BIBLIS" > "${RUN_DIR}/BIBLIS.stdout" 2>&1
  fi

  if [[ "${CASES}" == *"PNM"* ]]; then
    cp "${ROOT_DIR}/smpl/static/PNM" "${RUN_DIR}/PNM"
    "${KOMODO_BIN}" "${RUN_DIR}/PNM" > "${RUN_DIR}/PNM.stdout" 2>&1
  fi

  parse_metrics
}

parse_metrics() {
  local iaea3ds_out="${RUN_DIR}/IAEA3Ds.out"
  local a1t_coarse_out="${RUN_DIR}/A1t_coarse.out"
  local a1_out="${RUN_DIR}/A1.out"
  local a2_out="${RUN_DIR}/A2.out"
  local adjoint_out="${RUN_DIR}/adjoint.out"
  local biblis_out="${RUN_DIR}/BIBLIS.out"
  local pnm_out="${RUN_DIR}/PNM.out"

  {
    echo -e "metric\tvalue"

    if [[ -f "${iaea3ds_out}" ]]; then
      echo -e "iaea3ds_keff\t$(awk '/MULTIPLICATION EFFECTIVE \(K-EFF\)/{print $NF}' "${iaea3ds_out}" | tail -1)"
    fi

    if [[ -f "${a1_out}" ]]; then
      local a1_last_row
      a1_last_row="$(awk '
        BEGIN{flag=0}
        /Itr  Boron Conc\./{flag=1;next}
        flag && $1 ~ /^[0-9]+$/ {line=$0}
        flag && /Radial Power Distribution/{flag=0}
        END{print line}
      ' "${a1_out}")"
      echo -e "a1_final_boron\t$(awk '{print $2}' <<<"${a1_last_row}")"
      echo -e "a1_final_keff\t$(awk '{print $3}' <<<"${a1_last_row}")"
    fi

    if [[ -f "${a1t_coarse_out}" ]]; then
      local a1t_coarse_peak_row a1t_coarse_final_row
      local a1t_coarse_avg_fuel_temp a1t_coarse_max_fuel_temp
      local a1t_coarse_avg_mod_temp a1t_coarse_max_mod_temp
      local a1t_coarse_outlet_mod_temp a1t_coarse_avg_mod_dens a1t_coarse_outlet_mod_dens
      a1t_coarse_peak_row="$(awk '
        $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+\.[0-9]+$/ && $3 ~ /^-?[0-9]+\.[0-9]+$/ && $4 ~ /^[0-9]\.[0-9]+E[+-][0-9]+$/ {
          if (!seen || ($4+0) > max) { max=$4+0; row=$0; seen=1 }
        }
        END{print row}
      ' "${a1t_coarse_out}")"
      a1t_coarse_final_row="$(awk '$1 ~ /^[0-9]+$/ {row=$0} END{print row}' "${a1t_coarse_out}")"
      a1t_coarse_avg_fuel_temp="$(awk 'match($0, /:[[:space:]]*([0-9.]+)/, m) {print m[1]; exit}' <<<"$(awk '/^  AVERAGE FUEL TEMPERATURE/{print; exit}' "${a1t_coarse_out}")")"
      a1t_coarse_max_fuel_temp="$(awk 'match($0, /:[[:space:]]*([0-9.]+)/, m) {print m[1]; exit}' <<<"$(awk '/^  MAX FUEL CENTERLINE TEMPERATURE/{print; exit}' "${a1t_coarse_out}")")"
      a1t_coarse_avg_mod_temp="$(awk 'match($0, /:[[:space:]]*([0-9.]+)/, m) {print m[1]; exit}' <<<"$(awk '/^  AVERAGE MODERATOR TEMPERATURE/{print; exit}' "${a1t_coarse_out}")")"
      a1t_coarse_max_mod_temp="$(awk 'match($0, /:[[:space:]]*([0-9.]+)/, m) {print m[1]; exit}' <<<"$(awk '/^  MAXIMUM MODERATOR TEMPERATURE/{print; exit}' "${a1t_coarse_out}")")"
      a1t_coarse_outlet_mod_temp="$(awk 'match($0, /:[[:space:]]*([0-9.]+)/, m) {print m[1]; exit}' <<<"$(awk '/^  OUTLET MODERATOR TEMPERATURE/{print; exit}' "${a1t_coarse_out}")")"
      a1t_coarse_avg_mod_dens="$(awk 'match($0, /:[[:space:]]*([0-9.]+)/, m) {print m[1]; exit}' <<<"$(awk '/^  AVERAGE MODERATOR DENSITY/{print; exit}' "${a1t_coarse_out}")")"
      a1t_coarse_outlet_mod_dens="$(awk 'match($0, /:[[:space:]]*([0-9.]+)/, m) {print m[1]; exit}' <<<"$(awk '/^  OUTLET MODERATOR DENSITY/{print; exit}' "${a1t_coarse_out}")")"
      echo -e "a1t_coarse_peak_step\t$(awk '{print $1}' <<<"${a1t_coarse_peak_row}")"
      echo -e "a1t_coarse_peak_rel_power\t$(awk '{print $4}' <<<"${a1t_coarse_peak_row}")"
      echo -e "a1t_coarse_final_step\t$(awk '{print $1}' <<<"${a1t_coarse_final_row}")"
      echo -e "a1t_coarse_final_time\t$(awk '{print $2}' <<<"${a1t_coarse_final_row}")"
      echo -e "a1t_coarse_final_reactivity\t$(awk '{print $3}' <<<"${a1t_coarse_final_row}")"
      echo -e "a1t_coarse_final_rel_power\t$(awk '{print $4}' <<<"${a1t_coarse_final_row}")"
      echo -e "a1t_coarse_avg_fuel_temp\t${a1t_coarse_avg_fuel_temp}"
      echo -e "a1t_coarse_max_fuel_temp\t${a1t_coarse_max_fuel_temp}"
      echo -e "a1t_coarse_avg_mod_temp\t${a1t_coarse_avg_mod_temp}"
      echo -e "a1t_coarse_max_mod_temp\t${a1t_coarse_max_mod_temp}"
      echo -e "a1t_coarse_outlet_mod_temp\t${a1t_coarse_outlet_mod_temp}"
      echo -e "a1t_coarse_avg_mod_dens\t${a1t_coarse_avg_mod_dens}"
      echo -e "a1t_coarse_outlet_mod_dens\t${a1t_coarse_outlet_mod_dens}"
    fi

    if [[ -f "${a2_out}" ]]; then
      local a2_last_row
      a2_last_row="$(awk '
        BEGIN{flag=0}
        /Itr  Boron Conc\./{flag=1;next}
        flag && $1 ~ /^[0-9]+$/ {line=$0}
        flag && /Radial Power Distribution/{flag=0}
        END{print line}
      ' "${a2_out}")"
      echo -e "a2_final_boron\t$(awk '{print $2}' <<<"${a2_last_row}")"
      echo -e "a2_final_keff\t$(awk '{print $3}' <<<"${a2_last_row}")"
    fi

    if [[ -f "${adjoint_out}" ]]; then
      echo -e "adjoint_keff\t$(awk '/MULTIPLICATION EFFECTIVE \(K-EFF\)/{print $NF}' "${adjoint_out}" | tail -1)"
    fi

    if [[ -f "${biblis_out}" ]]; then
      echo -e "biblis_keff\t$(awk '/MULTIPLICATION EFFECTIVE \(K-EFF\)/{print $NF}' "${biblis_out}" | tail -1)"
    fi

    if [[ -f "${pnm_out}" ]]; then
      echo -e "pnm_keff\t$(awk '/MULTIPLICATION EFFECTIVE \(K-EFF\)/{print $NF}' "${pnm_out}" | tail -1)"
    fi
  } > "${METRICS_FILE}"

  cat <<EOF
Stage 0 run complete.
Artifacts: ${ARTIFACT_DIR}
Metrics:   ${METRICS_FILE}
EOF
}

get_metric() {
  local key="$1"
  awk -F'\t' -v k="${key}" '$1==k {print $2}' "${METRICS_FILE}"
}

assert_equal() {
  local name="$1"
  local actual="$2"
  local expected="$3"
  if [[ "${actual}" != "${expected}" ]]; then
    echo "FAIL: ${name}: actual=${actual}, expected=${expected}" >&2
    return 1
  fi
  echo "PASS: ${name}: ${actual}"
}

assert_close() {
  local name="$1"
  local actual="$2"
  local expected="$3"
  local tol="$4"
  local ok
  ok="$(awk -v a="${actual}" -v e="${expected}" -v t="${tol}" 'BEGIN {d=a-e; if (d<0) d=-d; print (d<=t) ? "1" : "0"}')"
  if [[ "${ok}" != "1" ]]; then
    echo "FAIL: ${name}: actual=${actual}, expected=${expected}, tol=${tol}" >&2
    return 1
  fi
  echo "PASS: ${name}: actual=${actual}, expected=${expected}, tol=${tol}"
}

check_stage0() {
  run_stage0

  local iaea3ds_keff a1_boron a1_keff
  local a1t_coarse_peak_step a1t_coarse_peak_relp a1t_coarse_final_step a1t_coarse_final_time a1t_coarse_final_react a1t_coarse_final_relp
  local a1t_coarse_avg_fuel_temp a1t_coarse_max_fuel_temp a1t_coarse_avg_mod_temp a1t_coarse_max_mod_temp
  local a1t_coarse_outlet_mod_temp a1t_coarse_avg_mod_dens a1t_coarse_outlet_mod_dens
  local a2_boron a2_keff
  local adjoint_keff biblis_keff pnm_keff

  iaea3ds_keff="$(get_metric iaea3ds_keff)"
  a1_boron="$(get_metric a1_final_boron)"
  a1_keff="$(get_metric a1_final_keff)"
  a1t_coarse_peak_step="$(get_metric a1t_coarse_peak_step)"
  a1t_coarse_peak_relp="$(get_metric a1t_coarse_peak_rel_power)"
  a1t_coarse_final_step="$(get_metric a1t_coarse_final_step)"
  a1t_coarse_final_time="$(get_metric a1t_coarse_final_time)"
  a1t_coarse_final_react="$(get_metric a1t_coarse_final_reactivity)"
  a1t_coarse_final_relp="$(get_metric a1t_coarse_final_rel_power)"
  a1t_coarse_avg_fuel_temp="$(get_metric a1t_coarse_avg_fuel_temp)"
  a1t_coarse_max_fuel_temp="$(get_metric a1t_coarse_max_fuel_temp)"
  a1t_coarse_avg_mod_temp="$(get_metric a1t_coarse_avg_mod_temp)"
  a1t_coarse_max_mod_temp="$(get_metric a1t_coarse_max_mod_temp)"
  a1t_coarse_outlet_mod_temp="$(get_metric a1t_coarse_outlet_mod_temp)"
  a1t_coarse_avg_mod_dens="$(get_metric a1t_coarse_avg_mod_dens)"
  a1t_coarse_outlet_mod_dens="$(get_metric a1t_coarse_outlet_mod_dens)"
  a2_boron="$(get_metric a2_final_boron)"
  a2_keff="$(get_metric a2_final_keff)"
  adjoint_keff="$(get_metric adjoint_keff)"
  biblis_keff="$(get_metric biblis_keff)"
  pnm_keff="$(get_metric pnm_keff)"

  assert_close "IAEA3Ds keff" "${iaea3ds_keff}" "${expect_iaea3ds_keff}" "${tol_keff}"
  assert_close "A1 final boron" "${a1_boron}" "${expect_a1_boron}" "${tol_boron}"
  assert_close "A1 final keff" "${a1_keff}" "${expect_a1_keff}" "${tol_keff}"

  assert_equal "A1t_coarse peak step" "${a1t_coarse_peak_step}" "${expect_a1t_coarse_peak_step}"
  assert_close "A1t_coarse peak relative power" "${a1t_coarse_peak_relp}" "${expect_a1t_coarse_peak_rel_power}" "${tol_rel_power}"
  assert_equal "A1t_coarse final step" "${a1t_coarse_final_step}" "${expect_a1t_coarse_final_step}"
  assert_close "A1t_coarse final time" "${a1t_coarse_final_time}" "${expect_a1t_coarse_final_time}" "1e-9"
  assert_close "A1t_coarse final reactivity" "${a1t_coarse_final_react}" "${expect_a1t_coarse_final_react}" "${tol_reactivity}"
  assert_close "A1t_coarse final relative power" "${a1t_coarse_final_relp}" "${expect_a1t_coarse_final_rel_power}" "${tol_rel_power}"
  assert_close "A1t_coarse average fuel temperature" "${a1t_coarse_avg_fuel_temp}" "${expect_a1t_coarse_avg_fuel_temp}" "1e-9"
  assert_close "A1t_coarse max fuel temperature" "${a1t_coarse_max_fuel_temp}" "${expect_a1t_coarse_max_fuel_temp}" "1e-9"
  assert_close "A1t_coarse average moderator temperature" "${a1t_coarse_avg_mod_temp}" "${expect_a1t_coarse_avg_mod_temp}" "1e-9"
  assert_close "A1t_coarse max moderator temperature" "${a1t_coarse_max_mod_temp}" "${expect_a1t_coarse_max_mod_temp}" "1e-9"
  assert_close "A1t_coarse outlet moderator temperature" "${a1t_coarse_outlet_mod_temp}" "${expect_a1t_coarse_outlet_mod_temp}" "1e-9"
  assert_close "A1t_coarse average moderator density" "${a1t_coarse_avg_mod_dens}" "${expect_a1t_coarse_avg_mod_dens}" "1e-9"
  assert_close "A1t_coarse outlet moderator density" "${a1t_coarse_outlet_mod_dens}" "${expect_a1t_coarse_outlet_mod_dens}" "1e-9"

  assert_close "A2 final boron" "${a2_boron}" "${expect_a2_boron}" "${tol_boron}"
  assert_close "A2 final keff" "${a2_keff}" "${expect_a2_keff}" "${tol_keff}"

  assert_close "Adjoint keff" "${adjoint_keff}" "${expect_adjoint_keff}" "${tol_keff}"
  assert_close "BIBLIS keff" "${biblis_keff}" "${expect_biblis_keff}" "${tol_keff}"
  assert_close "PNM keff" "${pnm_keff}" "${expect_pnm_keff}" "${tol_keff}"

  echo "Stage 0 baseline check: PASS"
}

cmd="${1:-}"
case "${cmd}" in
  run)
    run_stage0
    ;;
  check)
    check_stage0
    ;;
  *)
    usage
    exit 1
    ;;
esac
