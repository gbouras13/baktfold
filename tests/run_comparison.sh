#!/usr/bin/env bash
# run_comparison.sh
# Covers ALL non-euk baktfold test cases from test_integration.py.
# Runs each case in both:
#   baktfold     (dev: pholdlib-refactored)
#   baktfold_ref (bioconda v0.3.0 reference)
# then compares key output files ignoring timestamps.
#
# Usage (from any directory):
#   bash tests/run_comparison.sh
#
# All outputs → /tmp/baktfold_compare/

set -uo pipefail

trap 'echo "[EXIT] line=$LINENO status=$? sig=${_last_sig:-none}" >> /tmp/baktfold_compare_exit.log' EXIT
trap '_last_sig=HUP'  HUP
trap '_last_sig=TERM' TERM
trap '_last_sig=INT'  INT

CONDA="${CONDA_EXE:-$(command -v conda)}"
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DB="${REPO}/tests/test_data/baktfold_db"
TD="${REPO}/tests/test_data"
CASE_OUT="${TD}/outputs"        # each case writes to CASE_OUT/name
BASE="/tmp/baktfold_compare"
T=1                             # 1 thread for determinism

DEV_ENV="baktfold"
REF_ENV="baktfold_ref"
DEV_OUT="${BASE}/dev_outputs"
REF_OUT="${BASE}/ref_outputs"

# stable path for convert-prokka output (used by run_prokka)
PROKKA_JSON="${BASE}/assembly_prokka.json"

mkdir -p "$BASE" "$DEV_OUT" "$REF_OUT" "$CASE_OUT"

PASS=0
FAIL=0
ERRORS=()

# ── low-level helpers ────────────────────────────────────────────────────────

_run_env() {
    # _run_env ENV LABEL LOG CMD...
    local env="$1" label="$2" log="$3"; shift 3
    local cmd="$*"
    printf "  [%s] %s starting...\n" "$(date '+%H:%M:%S')" "$label"
    if $CONDA run -n "$env" bash -c "cd '$REPO' && $cmd" > "$log" 2>&1; then
        printf "  [%s] %s DONE ✓\n" "$(date '+%H:%M:%S')" "$label"
        return 0
    else
        printf "  [%s] %s FAILED ✗ (see %s)\n" "$(date '+%H:%M:%S')" "$label" "$log"
        ERRORS+=("$label FAILED — see $log")
        FAIL=$((FAIL + 1))
        return 1
    fi
}

_run_env_file() {
    # _run_env_file ENV LABEL LOG SCRIPT_FILE
    local env="$1" label="$2" log="$3" script_file="$4"
    printf "  [%s] %s starting...\n" "$(date '+%H:%M:%S')" "$label"
    if $CONDA run -n "$env" bash "$script_file" > "$log" 2>&1; then
        printf "  [%s] %s DONE ✓\n" "$(date '+%H:%M:%S')" "$label"
        return 0
    else
        printf "  [%s] %s FAILED ✗ (see %s)\n" "$(date '+%H:%M:%S')" "$label" "$log"
        ERRORS+=("$label FAILED — see $log")
        FAIL=$((FAIL + 1))
        return 1
    fi
}

_compare() {
    local name="$1" snap_dev="$2" snap_ref="$3"
    printf "  Comparing outputs...\n"
    if python3 "${REPO}/tests/compare_outputs.py" "$snap_dev" "$snap_ref" --cpu; then
        printf "  ✓ MATCH\n"
        PASS=$((PASS + 1))
    else
        printf "  ✗ DIFFER\n"
        ERRORS+=("$name: outputs differ")
        FAIL=$((FAIL + 1))
    fi
}

# ── case runners ─────────────────────────────────────────────────────────────

# Standard case: both envs run identical command, compare output dirs.
# CMD should write output to CASE_OUT/NAME.
run_case() {
    local name="$1"; shift
    local baktfold_cmd="$*"
    local snap_dir="${CASE_OUT}/${name}"
    local snap_dev="${DEV_OUT}/${name}"
    local snap_ref="${REF_OUT}/${name}"
    local dev_log="${BASE}/${name}_dev.log"
    local ref_log="${BASE}/${name}_ref.log"

    printf "\n══ %-50s ══\n" "$name"

    rm -rf "$snap_dir"; mkdir -p "$snap_dir"
    if _run_env "$DEV_ENV" "${name}[dev]" "$dev_log" "$baktfold_cmd"; then
        rm -rf "$snap_dev"; cp -r "$snap_dir" "$snap_dev"
    else
        return
    fi

    rm -rf "$snap_dir"; mkdir -p "$snap_dir"
    if _run_env "$REF_ENV" "${name}[ref]" "$ref_log" "$baktfold_cmd"; then
        rm -rf "$snap_ref"; cp -r "$snap_dir" "$snap_ref"
    else
        return
    fi

    _compare "$name" "$snap_dev" "$snap_ref"
}

# Script-file case: run a pre-written bash script (for tricky quoting).
run_case_via_script() {
    local name="$1" script_file="$2"
    local snap_dir="${CASE_OUT}/${name}"
    local snap_dev="${DEV_OUT}/${name}"
    local snap_ref="${REF_OUT}/${name}"
    local dev_log="${BASE}/${name}_dev.log"
    local ref_log="${BASE}/${name}_ref.log"

    printf "\n══ %-50s ══\n" "$name"

    rm -rf "$snap_dir"; mkdir -p "$snap_dir"
    if _run_env_file "$DEV_ENV" "${name}[dev]" "$dev_log" "$script_file"; then
        rm -rf "$snap_dev"; cp -r "$snap_dir" "$snap_dev"
    else
        return
    fi

    rm -rf "$snap_dir"; mkdir -p "$snap_dir"
    if _run_env_file "$REF_ENV" "${name}[ref]" "$ref_log" "$script_file"; then
        rm -rf "$snap_ref"; cp -r "$snap_dir" "$snap_ref"
    else
        return
    fi

    _compare "$name" "$snap_dev" "$snap_ref"
}

# Compare case: dev and ref each use their own predictions-dir snap.
# BASE_CMD should write output to CASE_OUT/NAME.
run_compare_case() {
    local name="$1" dev_pred="$2" ref_pred="$3"; shift 3
    local base_cmd="$*"
    local snap_dir="${CASE_OUT}/${name}"
    local snap_dev="${DEV_OUT}/${name}"
    local snap_ref="${REF_OUT}/${name}"
    local dev_log="${BASE}/${name}_dev.log"
    local ref_log="${BASE}/${name}_ref.log"

    printf "\n══ %-50s ══\n" "$name"

    rm -rf "$snap_dir"; mkdir -p "$snap_dir"
    if _run_env "$DEV_ENV" "${name}[dev]" "$dev_log" \
            "$base_cmd --predictions-dir '$dev_pred'"; then
        rm -rf "$snap_dev"; cp -r "$snap_dir" "$snap_dev"
    else
        return
    fi

    rm -rf "$snap_dir"; mkdir -p "$snap_dir"
    if _run_env "$REF_ENV" "${name}[ref]" "$ref_log" \
            "$base_cmd --predictions-dir '$ref_pred'"; then
        rm -rf "$snap_ref"; cp -r "$snap_dir" "$snap_ref"
    else
        return
    fi

    _compare "$name" "$snap_dev" "$snap_ref"
}

# Exit-code-only case: just verify both envs exit 0 (no file comparison).
run_exitcode_case() {
    local name="$1"; shift
    local baktfold_cmd="$*"
    local dev_log="${BASE}/${name}_dev.log"
    local ref_log="${BASE}/${name}_ref.log"

    printf "\n══ %-50s ══\n" "$name"

    if ! _run_env "$DEV_ENV" "${name}[dev]" "$dev_log" "$baktfold_cmd"; then
        return
    fi
    if _run_env "$REF_ENV" "${name}[ref]" "$ref_log" "$baktfold_cmd"; then
        printf "  ✓ both exited 0 (exit-code-only test)\n"
        PASS=$((PASS + 1))
    fi
}

# Single-file case: output is one file (not a directory).
# After each run the file is copied into a snap subdir for comparison.
run_single_file_case() {
    local name="$1" out_file="$2"; shift 2
    local baktfold_cmd="$*"
    local snap_dev="${DEV_OUT}/${name}"
    local snap_ref="${REF_OUT}/${name}"
    local dev_log="${BASE}/${name}_dev.log"
    local ref_log="${BASE}/${name}_ref.log"

    printf "\n══ %-50s ══\n" "$name"

    mkdir -p "$snap_dev" "$snap_ref"
    local fname; fname="$(basename "$out_file")"

    rm -f "$out_file"
    if _run_env "$DEV_ENV" "${name}[dev]" "$dev_log" "$baktfold_cmd"; then
        cp "$out_file" "${snap_dev}/${fname}"
    else
        return
    fi

    rm -f "$out_file"
    if _run_env "$REF_ENV" "${name}[ref]" "$ref_log" "$baktfold_cmd"; then
        cp "$out_file" "${snap_ref}/${fname}"
    else
        return
    fi

    _compare "$name" "$snap_dev" "$snap_ref"
}

# ═══════════════════════════════════════════════════════════════════════════
# TEST CASES
# ═══════════════════════════════════════════════════════════════════════════

# ── predict (3Di inference only, no foldseek) ─────────────────────────────

run_case "predict_json" \
    "baktfold predict \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/predict_json \
       -t $T -d $DB -f --cpu"

run_case "predict_embeddings" \
    "baktfold predict \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/predict_embeddings \
       -t $T -d $DB -f \
       --save-per-residue-embeddings --save-per-protein-embeddings --cpu"

# ── proteins-predict ──────────────────────────────────────────────────────

run_case "proteins_predict" \
    "baktfold proteins-predict \
       -i ${TD}/assembly.hypotheticals.faa \
       -o ${CASE_OUT}/proteins_predict \
       -t $T -d $DB -f --cpu"

run_case "proteins_predict_bakta_proteins" \
    "baktfold proteins-predict \
       -i ${TD}/assembly_bakta_proteins_output/assembly.hypotheticals.faa \
       -o ${CASE_OUT}/proteins_predict_bakta_proteins \
       -t $T -d $DB -f --cpu"

# ── run (3Di inference + foldseek + bakta annotation output) ─────────────

run_case "run" \
    "baktfold run \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/run \
       -t $T -d $DB -f --cpu"

run_case "run_no_fs_hits" \
    "baktfold run \
       -i ${TD}/SAMEA111266571.bakta.json \
       -o ${CASE_OUT}/run_no_fs_hits \
       -t $T -d $DB -e 1e-50 -f --cpu"

run_case "run_all" \
    "baktfold run \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/run_all \
       -t $T -d $DB -f -a --cpu"

run_case "run_fast" \
    "baktfold run \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/run_fast \
       -t $T -d $DB -f --fast --cpu"

# extra-foldseek-params: needs "--cov-mode 2" with embedded double quotes,
# so use a temp script to avoid bash -c quoting issues.
_extra_script="/tmp/baktfold_extra_foldseek_$$.sh"
cat > "$_extra_script" << ENDHEREDOC
#!/bin/bash
set -eo pipefail
cd '${REPO}'
exec baktfold run -i '${TD}/assembly_bakta_output/assembly.json' -o '${CASE_OUT}/run_extra_foldseek_params' -t ${T} -d '${DB}' -f --extra-foldseek-params "--cov-mode 2" --cpu
ENDHEREDOC
run_case_via_script "run_extra_foldseek_params" "$_extra_script"
rm -f "$_extra_script"

run_case "run_custom_db" \
    "baktfold run \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/run_custom_db \
       -t $T -d $DB \
       --custom-db ${TD}/custom_db/dummy_custom_db \
       -f --cpu"

run_case "run_custom_db_custom_annotations" \
    "baktfold run \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/run_custom_db_custom_annotations \
       -t $T -d $DB \
       --custom-db ${TD}/custom_db/dummy_custom_db \
       --custom-annotations ${TD}/custom_db/dummy_custom_db_annotations.tsv \
       -f --cpu"

# ── proteins (fasta → 3Di + foldseek, no bakta annotation output) ─────────

run_case "proteins" \
    "baktfold proteins \
       -i ${TD}/assembly.hypotheticals.faa \
       -o ${CASE_OUT}/proteins \
       -t $T -d $DB -f --cpu"

run_case "proteins_json" \
    "baktfold proteins \
       -i ${TD}/assembly_bakta_proteins_output_all/assembly.json \
       -o ${CASE_OUT}/proteins_json \
       -t $T -d $DB -f --cpu"

run_case "proteins_pipe" \
    "baktfold proteins \
       -i ${TD}/pipe.faa \
       -o ${CASE_OUT}/proteins_pipe \
       -t $T -d $DB -f --cpu"

# ── structure-based compare (no ProstT5 prerequisite) ────────────────────

run_case "compare_pdb" \
    "baktfold compare \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/compare_pdb \
       -t $T -d $DB \
       --structure-dir ${TD}/pdbs -f"

run_case "compare_cif" \
    "baktfold compare \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/compare_cif \
       -t $T -d $DB \
       --structure-dir ${TD}/cifs -f"

run_case "proteins_compare_pdb" \
    "baktfold proteins-compare \
       -i ${TD}/assembly.hypotheticals.faa \
       -o ${CASE_OUT}/proteins_compare_pdb \
       -t $T -d $DB \
       --structure-dir ${TD}/pdbs -f"

run_case "proteins_compare_cif" \
    "baktfold proteins-compare \
       -i ${TD}/assembly.hypotheticals.faa \
       -o ${CASE_OUT}/proteins_compare_cif \
       -t $T -d $DB \
       --structure-dir ${TD}/cifs -f"

# ── convert-prokka (single JSON file output) ──────────────────────────────

run_single_file_case "convert_prokka" "$PROKKA_JSON" \
    "baktfold convert-prokka \
       -i ${TD}/assembly_prokka_output/PROKKA_02192026.gbk \
       -o $PROKKA_JSON -f"

# ═══════════════════════════════════════════════════════════════════════════
# DEPENDENT CASES (require prior snap dirs from above)
# ═══════════════════════════════════════════════════════════════════════════

# compare: uses predict_json snap for --predictions-dir
run_compare_case "compare" \
    "${DEV_OUT}/predict_json" "${REF_OUT}/predict_json" \
    "baktfold compare \
       -i ${TD}/assembly_bakta_output/assembly.json \
       -o ${CASE_OUT}/compare \
       -t $T -d $DB -f"

# proteins_compare: uses proteins_predict snap for --predictions-dir
run_compare_case "proteins_compare" \
    "${DEV_OUT}/proteins_predict" "${REF_OUT}/proteins_predict" \
    "baktfold proteins-compare \
       -i ${TD}/assembly.hypotheticals.faa \
       -o ${CASE_OUT}/proteins_compare \
       -t $T -d $DB -f"

# run_prokka: uses prokka JSON from convert_prokka (deterministic, same for both)
run_case "run_prokka" \
    "baktfold run \
       -i $PROKKA_JSON \
       -o ${CASE_OUT}/run_prokka \
       -t $T -d $DB -f --cpu"

# ── summary ───────────────────────────────────────────────────────────────

echo ""
echo "══════════════════════════════════════════════════════"
printf "  Results: %d passed, %d failed\n" "$PASS" "$FAIL"
echo "══════════════════════════════════════════════════════"
if [ "${#ERRORS[@]}" -gt 0 ]; then
    echo "  Failures:"
    for e in "${ERRORS[@]}"; do
        echo "    ✗ $e"
    done
    exit 1
else
    echo "  All cases match ✓"
    exit 0
fi
