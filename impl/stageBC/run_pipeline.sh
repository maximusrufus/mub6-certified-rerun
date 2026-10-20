#!/bin/bash
# Stage B + Stage C pipeline of SPEC section 3 / 5 (execution.tex invocation order), starting from Stage A records.
#   run_pipeline.sh WORKDIR coarse0.bin [coarse1.bin ...]
# Environment overrides (defaults = the paper's parameters):
#   BIN=./build          directory with brun / crun
#   TAU_MID=.07 W_MID=.013   TAU_FINE=.016 W_FINE=.005   TAU_FB=.011 W_FB=.0038
#   CAP=5000000          graph step cap of the initial C pass          CAP_FINAL=5000000   cap of the final C pass
#   CSHARDS=4            shards of the final C pass
# Invariant: the run ends with "0 unresolved" (exit 0) or a residual list (exit 1; files RESIDUAL_*.bin in WORKDIR).
# Nothing is dropped silently: every B failure goes to a bad file, every C failure to a bad file.
set -u
W=$1; shift
BIN=${BIN:-./build}
TAU_MID=${TAU_MID:-.07}; W_MID=${W_MID:-.013}
TAU_FINE=${TAU_FINE:-.016}; W_FINE=${W_FINE:-.005}
TAU_FB=${TAU_FB:-.011}; W_FB=${W_FB:-.0038}
CAP=${CAP:-5000000}; CAP_FINAL=${CAP_FINAL:-5000000}; CSHARDS=${CSHARDS:-4}
export LC_ALL=C
mkdir -p "$W/batches"
rm -f "$W"/batches/*.bin "$W"/batches/fdone.txt "$W"/RESIDUAL_*.bin
RES=0

echo "== Stage B mid (tau=$TAU_MID w=$W_MID)"
"$BIN/brun" $TAU_MID $W_MID "$W/mid.bin" "$W/midbad.bin" - "$@" || { echo "B mid: failed parents (residual)"; cp "$W/midbad.bin" "$W/RESIDUAL_Bmid.bin"; RES=1; }

echo "== Stage B fine (tau=$TAU_FINE w=$W_FINE), batches"
"$BIN/brun" $TAU_FINE $W_FINE "$W/fine.bin" "$W/finebad.bin" "$W/batches/" "$W/mid.bin" || { echo "B fine: failed parents (residual)"; cp "$W/finebad.bin" "$W/RESIDUAL_Bfine.bin"; RES=1; }
[ -f "$W/batches/fdone.txt" ] || { echo "no fdone.txt marker"; exit 2; }
read NB NT < "$W/batches/fdone.txt"
echo "   $NB batches, $NT roots (marker)"

echo "== Stage C initial (one process per batch)"
BADFILES=()
for ((j = 0; j < NB; j++)); do
  "$BIN/crun" 12 4 1 .042 1 1 $CAP 0 1 "$W/batches/bad$j.bin" "$W/batches/f$j.bin" > /dev/null 2> "$W/batches/c$j.log"
  rc=$?
  if [ $rc -ge 2 ]; then echo "crun batch $j: error (exit $rc)"; exit 2; fi
  SZ=$(stat -c %s "$W/batches/bad$j.bin")
  if [ "$SZ" -gt 8 ]; then BADFILES+=("$W/batches/bad$j.bin"); fi
done
echo "   batches with unresolved records: ${#BADFILES[@]}"

if [ ${#BADFILES[@]} -gt 0 ]; then
  echo "== Stage B fallback (tau=$TAU_FB w=$W_FB) on the C failures (C-locale filename order)"
  SORTED=($(printf '%s\n' "${BADFILES[@]}" | sort))
  "$BIN/brun" $TAU_FB $W_FB "$W/fb1.bin" "$W/fbbad.bin" - "${SORTED[@]}" || { echo "B fallback: failed parents (residual)"; cp "$W/fbbad.bin" "$W/RESIDUAL_Bfb.bin"; RES=1; }
  echo "== Stage C final ($CSHARDS shards)"
  for ((i = 0; i < CSHARDS; i++)); do
    "$BIN/crun" 12 4 1 .042 1 1 $CAP_FINAL $i $CSHARDS "$W/finalbad$i.bin" "$W/fb1.bin" > /dev/null 2> "$W/final$i.log"
    rc=$?
    if [ $rc -ge 2 ]; then echo "final crun shard $i: error (exit $rc)"; exit 2; fi
    SZ=$(stat -c %s "$W/finalbad$i.bin")
    if [ "$SZ" -gt 8 ]; then echo "   shard $i: unresolved records remain"; cp "$W/finalbad$i.bin" "$W/RESIDUAL_C$i.bin"; RES=1; fi
  done
fi
if [ $RES -eq 0 ]; then echo "RESULT: 0 unresolved"; else echo "RESULT: residual list in $W/RESIDUAL_*.bin"; fi
exit $RES
