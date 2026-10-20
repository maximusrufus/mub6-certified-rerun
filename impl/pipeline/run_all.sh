#!/bin/bash
# Official pipeline driver:  seed (regenerated) -> Stage A (arun) -> S5/S6 verification -> Stage B (brun) -> Stage C (crun)
#                            -> fallback loop to a fixed point -> MANIFEST.
#   run_all.sh [WORKDIR=pipeline/work]
# Environment: STRIDE=k keep every k-th seed line (1 = all; 100 = 1% smoke sample)    THREADS=16
#   SAMPLE_MOD=100  S6 samples 1/SAMPLE_MOD of all EXCLUDED Stage A cells (splitmix64 hash, deterministic)
#   BC_BUDGET=0     seconds allowed for stage B+C (0 = unlimited). On expiry the B/C process group is killed and the run
#                   ends as INCOMPLETE (a residual with a reason), never as "0 unresolved".
#   MANIFEST=pipeline/MANIFEST.txt
# Last stdout line: exactly "0 unresolved" (exit 0), or "unresolved (residual list above)" after the reasons (exit 1).
# Exit 2 = tool/format/verification failure (self-test, hand-off format, S5 re-certification, S6 disagreement).
set -u
export LC_ALL=C
P=$(cd "$(dirname "$0")" && pwd); I=$(dirname "$P")
W=${1:-$P/work}; STRIDE=${STRIDE:-1}; THREADS=${THREADS:-16}; SAMPLE_MOD=${SAMPLE_MOD:-100}; BC_BUDGET=${BC_BUDGET:-0}
MANIFEST=${MANIFEST:-$P/MANIFEST.txt}
BIN=$P/bin; A=$I/stageA; BC=$I/stageBC
TAU_FINE=${TAU_FINE:-.016}; TAU_FB=${TAU_FB:-.011}; W_FB=${W_FB:-.0038}; CAP_FINAL=${CAP_FINAL:-5000000}; CSHARDS=${CSHARDS:-4}; MAXROUNDS=${MAXROUNDS:-20}
mkdir -p "$BIN" "$W"; W=$(cd "$W" && pwd)
T0=$(date +%s.%N); TL=$T0
tick(){ N=$(date +%s.%N); printf "[time] %s: %.1f s (cumulative %.1f s)\n" "$1" "$(echo "$N - $TL" | bc)" "$(echo "$N - $T0" | bc)"; TL=$N; }
# S3: every official binary is built with -DIV_CHECK_ROUNDING (per-operation FE_UPWARD assertion).  Never -ffast-math.
CXXF="-O2 -frounding-math -ffp-contract=off -fopenmp -Wall -Wextra -DIV_CHECK_ROUNDING"
need(){ out=$1; shift; [ -x "$out" ] || return 0; for f in "$@"; do [ "$f" -nt "$out" ] && return 0; done; return 1; }
CORE="$I/core/interval.hpp $I/core/mub_ops.hpp $I/core/selftest.hpp $I/core/tensor_c.hpp"
if need "$BIN/seed" $I/seed/*.cpp $I/seed/*.hpp $CORE; then (cd $I/seed && g++ -std=gnu++17 $CXXF -o "$BIN/seed" seed.cpp) || { echo "build seed failed"; exit 2; }; fi
if need "$BIN/arun" $A/*.cpp $A/*.hpp $CORE; then (cd $A && g++ -std=c++20 $CXXF -o "$BIN/arun" arun.cpp seedproc.cpp root.cpp cover.cpp atlas.cpp) || { echo "build arun failed"; exit 2; }; fi
if need "$BIN/test_stagea" $A/*.cpp $A/*.hpp $A/tests/test_stagea.cpp $CORE; then (cd $A && g++ -std=c++20 $CXXF -o "$BIN/test_stagea" tests/test_stagea.cpp seedproc.cpp root.cpp cover.cpp atlas.cpp) || { echo "build test_stagea failed"; exit 2; }; fi
for t in brun crun; do
  # stageBC/build.sh already passes -DIV_CHECK_ROUNDING
  if need "$BIN/$t" $BC/*.cpp $BC/*.hpp $CORE; then bash $BC/build.sh $t.cpp "$BIN/$t" || { echo "build $t failed"; exit 2; }; fi
done
tick "build"
rm -rf "$W"/* 2>/dev/null; mkdir -p "$W"

# ---- S1: regenerate the seed list from the hashed seed binary; record provenance
echo "== seed (regenerated, interval build)"
"$BIN/seed" --out "$W/seeds_gen.txt" --quiet > "$W/seed.log" 2>&1 || { echo "seed failed (self-test or run)"; tail -3 "$W/seed.log"; exit 2; }
grep -h "selftest" "$W/seed.log" | head -1
NGEN=$(wc -l < "$W/seeds_gen.txt"); GENSHA=$(sha256sum "$W/seeds_gen.txt" | cut -d' ' -f1)
if cmp -s "$W/seeds_gen.txt" "$I/seed/seeds.txt"; then SAME=yes; else SAME=NO; fi
echo "seed provenance: $NGEN centres, sha256 $GENSHA, identical to committed seed/seeds.txt: $SAME"
[ "$SAME" = yes ] || echo "WARNING: regenerated seeds differ from seed/seeds.txt (the regenerated list is used)"
tick "seed regeneration"
SEEDS=$W/seeds_gen.txt
if [ "$STRIDE" -gt 1 ]; then awk -v k="$STRIDE" '(NR-1)%k==0' "$SEEDS" > "$W/seeds_in.txt"; else cp "$SEEDS" "$W/seeds_in.txt"; fi
NS=$(wc -l < "$W/seeds_in.txt"); echo "== seeds used: $NS lines (stride $STRIDE)"
awk 'NF!=10 || $1<0 || $1>4 {bad++} {for(i=2;i<=10;i++) if($i<0||$i>111||$i!~/^[0-9]+$/) bad++} END{exit bad>0}' "$W/seeds_in.txt" \
  || { echo "seed format incompatible with arun (expect 's b0..b8', b<112)"; exit 2; }

# ---- Stage A
echo "== Stage A (IV_CHECK_ROUNDING build, S6 sampling 1/$SAMPLE_MOD of excluded cells)"
"$BIN/arun" "$W/seeds_in.txt" 0 1 "$W/A.bin" "$W/Abad.txt" --threads $THREADS --sample-excl $SAMPLE_MOD "$W/excl_sample.txt" > "$W/A.log" 2>&1; rc=$?
grep -h "selftest" "$W/A.log" | head -1
[ $rc -ge 2 ] && { echo "arun error (exit $rc)"; tail -5 "$W/A.log"; exit 2; }
NABAD=$(wc -l < "$W/Abad.txt"); tail -n 1 "$W/A.log" | cut -c1-200; tick "stage A"

# ---- S5: re-certify ALL Stage A records (C++ checkRoot) + flint sample; hand-off format check
echo "== S5 re-certification of Stage A records"
"$BIN/test_stagea" --verify "$W/A.bin" > "$W/S5_checkRoot.log" 2>&1 || { tail -5 "$W/S5_checkRoot.log"; echo "S5: checkRoot re-certification FAILED"; exit 2; }
tail -2 "$W/S5_checkRoot.log" | tr '\n' ' '; echo
python3 "$A/tools/verify_records.py" "$W/A.bin" 40 > "$W/S5_flint.log" 2>&1 || { tail -5 "$W/S5_flint.log"; echo "S5: python-flint record check FAILED"; exit 2; }
tail -1 "$W/S5_flint.log"
python3 "$P/handoff_check.py" "$W/A.bin" || { echo "A->B hand-off format check FAILED"; exit 2; }
tick "S5 + hand-off"

# ---- S6: python-flint re-verification of the sampled excluded cells
echo "== S6 excluded-cell sample re-verification (python-flint)"
python3 "$P/verify_excluded.py" "$W/seeds_in.txt" "$W/excl_sample.txt" $THREADS > "$W/S6.log" 2>&1; rc=$?
tail -3 "$W/S6.log" | sed 's/^/S6: /'
[ $rc -ne 0 ] && { echo "S6 FAILED (disagreement or empty sample)"; exit 2; }
tick "S6"

# ---- Stage B + C (stageBC/run_pipeline.sh unmodified: mid -> fine+batches -> C per batch -> B fallback -> final C)
echo "== Stage B + C"
export BIN OMP_NUM_THREADS=$THREADS
BCRC=0
if [ "$BC_BUDGET" -gt 0 ]; then
  setsid bash $BC/run_pipeline.sh "$W" "$W/A.bin" > "$W/BC.log" 2>&1 &
  BCP=$!
  el=0
  while kill -0 $BCP 2>/dev/null; do
    sleep 5; el=$((el+5))
    if [ $el -ge $BC_BUDGET ]; then kill -TERM -- -$BCP 2>/dev/null; sleep 2; pkill -KILL -g $BCP 2>/dev/null; BCRC=124; break; fi
  done
  if [ $BCRC -ne 124 ]; then wait $BCP; BCRC=$?; fi
else
  bash $BC/run_pipeline.sh "$W" "$W/A.bin" > "$W/BC.log" 2>&1; BCRC=$?
fi
grep -E "^(==|   |B |RESULT|crun|final|no fdone)" "$W/BC.log" | head -30
LASTB=$(grep -E "^\[B\]" "$W/BC.log" 2>/dev/null | tail -1)
[ $BCRC -eq 2 ] && { echo "stage B/C error (exit 2)"; exit 2; }
tick "stage B/C round 1 (rc=$BCRC)"

# ---- fallback loop to a fixed point on the final-C residual (same fallback parameters; fixed point = residual unchanged)
round=1; prev=""
while [ $BCRC -ne 124 ] && ls "$W"/RESIDUAL_C*.bin > /dev/null 2>&1 && [ $round -lt $MAXROUNDS ]; do
  cur=$(cat "$W"/RESIDUAL_C*.bin | sha256sum | cut -d' ' -f1)
  if [ "$cur" = "$prev" ]; then echo "fallback fixed point reached (residual unchanged after round $round)"; break; fi
  prev=$cur; round=$((round+1)); echo "== fallback round $round"
  SORTED=($(ls "$W"/RESIDUAL_C*.bin | sort))
  "$BIN/brun" $TAU_FB $W_FB "$W/fbr$round.bin" "$W/fbrbad$round.bin" - "${SORTED[@]}" 2> "$W/fbr$round.log" || { echo "B fallback r$round: failed parents"; cp "$W/fbrbad$round.bin" "$W/RESIDUAL_Bfb$round.bin"; }
  rm -f "$W"/RESIDUAL_C*.bin
  for ((i = 0; i < CSHARDS; i++)); do
    "$BIN/crun" 12 4 1 .042 1 1 $CAP_FINAL $i $CSHARDS "$W/rbad${round}_$i.bin" "$W/fbr$round.bin" > /dev/null 2> "$W/rc${round}_$i.log"
    rc=$?; [ $rc -ge 2 ] && { echo "crun round $round shard $i: error"; exit 2; }
    [ "$(stat -c %s "$W/rbad${round}_$i.bin")" -gt 8 ] && cp "$W/rbad${round}_$i.bin" "$W/RESIDUAL_C${round}_$i.bin"
  done
done
tick "fallback loop ($round round(s))"

# ---- S5 (B/C part): structural pipeline check on a COMPLETED B/C work dir
if [ $BCRC -ne 124 ] && [ -f "$W/fine.bin" ]; then
  python3 $BC/tests/pipeline_check.py "$W" $TAU_FINE > "$W/S5_pipeline_check.log" 2>&1 || { tail -5 "$W/S5_pipeline_check.log"; echo "S5: pipeline_check FAILED"; exit 2; }
  echo "S5 pipeline_check: $(tail -1 "$W/S5_pipeline_check.log")"
fi

# ---- manifest
{ echo "# SHA-256 manifest, generated $(date -u +%FT%TZ); stride=$STRIDE threads=$THREADS sample_mod=$SAMPLE_MOD bc_budget=$BC_BUDGET"
  echo "# seed provenance: $NGEN centres regenerated by bin/seed; sha256 $GENSHA; identical to seed/seeds.txt: $SAME"
  echo "# S2 self-test lines logged: seed=$(grep -c 'startup selftest ok' "$W/seed.log") arun=$(grep -c 'startup selftest ok' "$W/A.log") brun/crun=$(cat "$W"/BC.log "$W"/batches/c*.log 2>/dev/null | grep -c 'startup selftest ok')"
  echo "# S6: $(tail -1 "$W/S6.log")   S5: $(tail -2 "$W/S5_checkRoot.log" | tr '\n' ' ')"
  for f in seed arun brun crun test_stagea; do echo "$(sha256sum "$BIN/$f" | cut -d' ' -f1)  bin/$f"; done
  echo "$(sha256sum "$I/seed/seeds.txt" | cut -d' ' -f1)  seed/seeds.txt (committed)"
  ( cd "$W" && for f in seeds_gen.txt seeds_in.txt Abad.txt A.bin excl_sample.txt mid.bin midbad.bin fine.bin finebad.bin fb1.bin fbbad.bin fbr*.bin finalbad*.bin rbad*.bin batches/f*.bin batches/bad*.bin RESIDUAL_*.bin; do
      [ -f "$f" ] && echo "$(sha256sum "$f" | cut -d' ' -f1)  out/$f"; done )
} > "$MANIFEST"
echo "manifest: $MANIFEST ($(grep -vc '^#' "$MANIFEST") hashes)"

# ---- verdict
RES=0
if [ "$NABAD" -gt 0 ]; then echo "residual: $NABAD seed(s) unresolved at Stage A (see $W/Abad.txt): $(head -3 "$W/Abad.txt" | tr '\n' '|')"; RES=1; fi
if [ $BCRC -eq 124 ]; then
  RES=1
  NREC=$(python3 -c "import struct;print(struct.unpack_from('<Q',open('$W/A.bin','rb').read(),0)[0])")
  echo "residual: stage B/C INCOMPLETE: not finished within BC_BUDGET=${BC_BUDGET}s (last: $LASTB); no B/C verdict exists for the $NREC Stage A records"
fi
for f in "$W"/RESIDUAL_*.bin; do
  [ -f "$f" ] || continue; RES=1; n=$(basename $f)
  case $n in RESIDUAL_Bmid*) why="Stage B mid failed parents";; RESIDUAL_Bfine*) why="Stage B fine failed parents";; RESIDUAL_Bfb*) why="Stage B fallback failed parents";; *) why="Stage C unresolved after fallback fixed point";; esac
  echo "residual: $n ($(stat -c %s $f) bytes): $why"
done
if [ $RES -eq 0 ]; then echo "0 unresolved"; exit 0; else echo "unresolved (residual list above)"; exit 1; fi
