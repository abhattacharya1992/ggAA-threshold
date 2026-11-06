#!/usr/bin/env bash
# Rebuild the threshold series of the two-loop Type II remainders from the published files.
#
#   bash reproduce_series.sh --inputs DIR --nf1-inputs DIR --work DIR [--package DIR]
#                            [--kernel PATH] [--block-seconds N]
#
#   --inputs         directory with Finite_2Loop_P1 and Finite_2Loop_P2 of arXiv:1909.08993
#   --nf1-inputs     directory with P1_nf1.wl and P1_nf1_compact.wl
#   --work           working directory; stages that finished are skipped when rerun
#   --package        ggAA-threshold package; the result is compared with its data/construction
#   --kernel         Wolfram kernel (default: math, WolframKernel or wolfram on PATH)
#   --block-seconds  time limit per exact rational operation (default 3600)
#
# Stages: nf2, p1_nf1_construction, p1_nf1_check, remaining, compare.
# Needs Mathematica 13.2, g++ (C++17) with GMP or Boost.Multiprecision headers,
# and Python 3.9+ with sympy and mpmath.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
inputs='' nf1='' work='' package='' kernel='' block_seconds=3600

usage() { sed -n '2,17p' "$0" | sed 's/^# \{0,1\}//'; }
die() { echo "ERROR: $*" >&2; exit 1; }
say() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

while (($#)); do
    case $1 in
        --inputs) inputs=${2:?}; shift 2 ;;
        --nf1-inputs) nf1=${2:?}; shift 2 ;;
        --work) work=${2:?}; shift 2 ;;
        --package) package=${2:?}; shift 2 ;;
        --kernel) kernel=${2:?}; shift 2 ;;
        --block-seconds) block_seconds=${2:?}; shift 2 ;;
        -h|--help) usage; exit 0 ;;
        *) usage; die "unknown option $1" ;;
    esac
done
[[ -n $inputs && -n $nf1 && -n $work ]] || { usage; exit 2; }

# ---------------------------------------------------------------- prerequisites
for c in g++ python3 sha256sum setsid; do
    command -v "$c" > /dev/null || die "$c not found"
done
if [[ -z $kernel ]]; then
    for c in math WolframKernel wolfram; do
        if command -v "$c" > /dev/null; then kernel=$(command -v "$c"); break; fi
    done
fi
[[ -n $kernel ]] || die "no Wolfram kernel found; use --kernel"
python3 -c 'import sys, sympy, mpmath; assert sys.version_info >= (3, 9)' \
    || die "Python 3.9+ with sympy and mpmath is required"

check_hash() {
    [[ -r $1 ]] || die "cannot read $1"
    [[ $(sha256sum "$1" | cut -d' ' -f1) == "$2" ]] || die "$1 differs from the expected file (SHA-256)"
}
inputs=$(cd "$inputs" && pwd)
nf1=$(cd "$nf1" && pwd)
check_hash "$inputs/Finite_2Loop_P1" 0ad929948dad92936c5bf8a7ae00ef90b8956d0e4412e8fc74a68c77a8c315b1
check_hash "$inputs/Finite_2Loop_P2" d7f76dbecc9c3aa4e41cf8c6d6b2d8191626c48d6e84a9e51e83ff5d342bd8e0
check_hash "$nf1/P1_nf1.wl" 0da89f4e7a22e433478b3e810b13cd720cf6c5eb29e168f771d1acd200be3815
check_hash "$nf1/P1_nf1_compact.wl" af362f3a8757a4335570347ccbbd781f87d060df87d163ab123bf37ac21958f4

mkdir -p "$work"
work=$(cd "$work" && pwd)
exec > >(tee -a "$work/reproduce.log") 2>&1
say "kernel: $kernel"
say "work directory: $work"

finished() { [[ -f $work/.done_$1 ]]; }
mark() { date -Is > "$work/.done_$1"; }
run_kernel() { "$kernel" -noinit -noprompt -script "$1" < /dev/null > "$2" 2>&1; }

# The exact rational code uses GMP if it is installed, Boost.Multiprecision otherwise.
cxx_exact=()
cxx_libs=()
if printf '#include <gmpxx.h>\nint main(){mpz_class a=3;return a*a==9?0:1;}\n' \
        | g++ -std=c++17 -x c++ - -lgmpxx -lgmp -o "$work/.gmp_probe" 2> /dev/null; then
    cxx_exact=(-DGGAA_USE_GMP)
    cxx_libs=(-lgmpxx -lgmp)
fi
rm -f "$work/.gmp_probe"

# ---------------------------------------------------------------- 1. Nf^2 coefficients
if ! finished nf2; then
    say "nf2: extracting and expanding the Nf^2 coefficients"
    d=$work/nf2
    rm -rf "$d"
    mkdir -p "$d/payload/source" "$d/payload/data" "$d/payload/checks" "$d/results" "$d/native"
    cp "$here"/nf2/* "$d/payload/source/"
    g++ -O2 -std=c++17 "$here/nf2/extract_nf.cpp" -o "$d/native/extract_nf"
    g++ -O2 -std=c++17 "$here/nf2/check_nf_extraction.cpp" -o "$d/native/check_nf_extraction"
    for i in 1 2; do
        "$d/native/extract_nf" "$inputs/Finite_2Loop_P$i" "$d/payload/data/P$i"
        "$d/native/check_nf_extraction" "$inputs/Finite_2Loop_P$i" "$d/payload/data/P${i}_nf2.wl"
    done
    cp "$here"/nf2/*.py "$d/payload/data/"
    (
        cd "$d/payload/data"
        python3 parse_nf2.py
        python3 collect_nf2_normalized.py P1
        python3 collect_nf2_normalized.py P2
        python3 derive_nf2.py
        python3 validate_nf2_python.py
    )
    mv "$d/payload/data/python_nf2_checks.json" "$d/payload/checks/"
    GGAA_NF2_RUN=$d GGAA_NF2_NO_PLOTS=1 run_kernel "$d/payload/source/driver.wl" "$d/kernel.log" \
        || die "Nf^2 driver failed; see $d/kernel.log"
    for p in P1 P2; do
        for o in 4 6; do
            [[ -s $d/results/${p}_series$o.wl ]] || die "Nf^2 driver did not write ${p}_series$o.wl"
        done
    done
    mark nf2
fi

# ---------------------------------------------------------------- 2. Nf coefficient of P1
if ! finished p1_nf1_construction; then
    say "p1_nf1_construction: threshold expansion of the Nf coefficient of P1"
    d=$work/p1_nf1_construction
    rm -rf "$d"
    mkdir -p "$d/input" "$d/source" "$d/results" "$d/checkpoints"
    cp "$nf1/P1_nf1.wl" "$nf1/P1_nf1_compact.wl" "$d/input/"
    cp "$here"/p1_nf1/construction/* "$d/source/"
    if ! GGAA_NF1_THRESHOLD_RUN=$d GGAA_NF1_NO_PLOTS=1 run_kernel "$d/source/driver.wl" "$d/kernel.log"; then
        # The original run also stopped here, in its own zero test, after writing both series.
        say "construction driver stopped after the expansion; the series are checked in the next stage"
    fi
    for o in 4 6; do
        [[ -s $d/results/P1_nf1_series$o.wl ]] || die "P1_nf1_series$o.wl was not written; see $d/kernel.log"
    done
    mark p1_nf1_construction
fi

if ! finished p1_nf1_check; then
    say "p1_nf1_check: exact and numerical checks of the P1 Nf series"
    d=$work/p1_nf1_check
    c=$work/p1_nf1_construction
    rm -rf "$d"
    mkdir -p "$d/input" "$d/source" "$d/results" "$d/checkpoints"
    cp "$nf1/P1_nf1.wl" "$nf1/P1_nf1_compact.wl" "$d/input/"
    for o in 4 6; do cp "$c/results/P1_nf1_series$o.wl" "$d/input/saved_P1_nf1_series$o.wl"; done
    cp "$here"/p1_nf1/check/* "$d/source/"
    GGAA_NF1_THRESHOLD_RUN=$d GGAA_NF1_REPAIR_SOURCE=$c GGAA_NF1_NO_PLOTS=1 \
        run_kernel "$d/source/driver.wl" "$d/kernel.log" || die "P1 Nf checks failed; see $d/kernel.log"
    for o in 4 6; do
        [[ -s $d/results/P1_nf1_series$o.wl ]] || die "check driver did not write P1_nf1_series$o.wl"
    done
    mark p1_nf1_check
fi

# ---------------------------------------------------------------- 3. remaining coefficients
d=$work/remaining
if ! finished remaining_setup; then
    say "remaining: collecting the source expressions into blocks"
    rm -rf "$d"
    mkdir -p "$d/input" "$d/known" "$d/collected" "$d/source" "$d/results" "$d/checkpoints"
    cp "$inputs/Finite_2Loop_P1" "$inputs/Finite_2Loop_P2" "$d/input/"
    for p in P1 P2; do
        for o in 4 6; do cp "$work/nf2/results/${p}_series$o.wl" "$d/known/${p}_nf2_p$o.wl"; done
    done
    for o in 4 6; do cp "$work/p1_nf1_check/results/P1_nf1_series$o.wl" "$d/known/P1_nf1_p$o.wl"; done
    cp "$here"/remaining/*.wl "$here"/remaining/*.py "$here"/remaining/*.cpp "$here"/remaining/*.hpp "$d/source/"
    cp "$here"/remaining/first_run/*.cpp "$d/source/"
    for name in collect_complete check_collection; do
        g++ -O2 -std=c++17 "$d/source/$name.cpp" -o "$d/source/$name"
    done
    for name in rational_dag check_fixed_reference jet_dag reference_dag; do
        g++ -O2 -std=c++17 ${cxx_exact[@]+"${cxx_exact[@]}"} "$d/source/$name.cpp" \
            ${cxx_libs[@]+"${cxx_libs[@]}"} -o "$d/source/$name"
    done
    (
        # the collector parses deeply nested expressions
        ulimit -s 65536 2> /dev/null || true
        for p in P1 P2; do
            "$d/source/collect_complete" "$d/input/Finite_2Loop_$p" "$d/collected/$p"
            "$d/source/check_collection" "$d/input/Finite_2Loop_$p" "$d/collected/$p"
        done
    )
    mark remaining_setup
fi

if ! finished remaining_rational; then
    say "remaining: exact rational series of all blocks (resumes where it stopped)"
    ( ulimit -s 65536 2> /dev/null || true
      python3 "$d/source/prepare_rational_series.py" "$d" --block-seconds "$block_seconds" )
    mark remaining_rational
fi

if ! finished remaining; then
    say "remaining: Nf^0 coefficients of P1 and P2, Nf coefficient of P2, assembly"
    GGAA_UNIFIED_RUN=$d GGAA_UNIFIED_PROJECT=$work GGAA_UNIFIED_NO_PLOTS=1 \
        setsid "$kernel" -noinit -noprompt -script "$d/source/driver.wl" < /dev/null > "$d/kernel.log" 2>&1 &
    pid=$!
    # Stage 09 of the driver is a long numerical validation that the original run did not
    # complete; the series are final once stage 08 has written the assembled form factors.
    while kill -0 "$pid" 2> /dev/null; do
        stage=$(head -n 1 "$d/current_stage.txt" 2> /dev/null || true)
        if [[ $stage == 09_full_numerical_validation && -s $d/results/P2_full_series.wl ]]; then
            say "series assembled; stopping before the driver's numerical validation"
            kill -TERM -- "-$pid" 2> /dev/null || kill -TERM "$pid" 2> /dev/null || true
            break
        fi
        sleep 30
    done
    wait "$pid" 2> /dev/null || true
    for f in P1_nf0 P2_nf0 P2_nf1; do
        for o in 4 6; do
            [[ -s $d/results/${f}_series$o.wl ]] || die "${f}_series$o.wl was not written; see $d/kernel.log"
        done
    done
    [[ -s $d/results/P1_full_series.wl && -s $d/results/P2_full_series.wl ]] \
        || die "assembly incomplete; see $d/kernel.log"
    mark remaining
fi

# ---------------------------------------------------------------- 4. comparison
if [[ -n $package ]]; then
    say "comparing with $package/data/construction"
    python3 "$here/compare_with_package.py" --work "$work" --package "$package"
fi
say "done; series in $work/remaining/results and $work/remaining/known"
