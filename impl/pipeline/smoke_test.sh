#!/bin/bash
# Smoke test: deterministic 1% sample (every 100th seed line) end to end.
# Exit 0 iff (a) run_all.sh ends in a clean status ("0 unresolved", or a residual list with reasons), and (b) the self-tests,
# S5 re-certification, S6 excluded-cell sample (0 disagreements) and the manifest all passed, and (c) the synthetic
# Stage-A-like B/C pipeline test (stageBC/tests/pipeline_test.sh) passed.
# SMOKE_BC_BUDGET (default 600 s) bounds the B/C stage on the REAL Stage A records (unbounded it takes days, see report);
# an expired budget is reported as an INCOMPLETE residual, never as success.  SMOKE_BC_BUDGET=0 = unbounded.
P=$(cd "$(dirname "$0")" && pwd); I=$(dirname "$P")
S=$(date +%s.%N)
STRIDE=100 BC_BUDGET=${SMOKE_BC_BUDGET:-600} bash "$P/run_all.sh" "$P/work_smoke" | tee "$P/work_smoke.log"
rc=${PIPESTATUS[0]}
echo "[time] run_all total: $(echo "$(date +%s.%N) - $S" | bc) s"
if [ $rc -ge 2 ]; then echo "SMOKE FAIL: run_all exit $rc (tool/format/verification failure)"; exit 1; fi
L=$P/work_smoke.log
if [ $rc -eq 0 ]; then
  tail -n 1 $L | grep -q '^0 unresolved$' || { echo "SMOKE FAIL: no 0-unresolved line"; exit 1; }
else
  grep -q '^residual: .*: .\+' $L || { echo "SMOKE FAIL: residual without reasons"; exit 1; }
fi
grep -q 'startup selftest ok' $L || { echo "SMOKE FAIL: S2 self-test line missing"; exit 1; }
grep -q '^S6: disagreements: 0$' $L || { echo "SMOKE FAIL: S6 did not report 0 disagreements"; exit 1; }
grep -q 'ALL OK' $P/work_smoke/S5_checkRoot.log || { echo "SMOKE FAIL: S5 checkRoot"; exit 1; }
[ "$(grep -c '^[0-9a-f]\{64\}  ' $P/MANIFEST.txt)" -ge 8 ] || { echo "SMOKE FAIL: manifest missing hashes"; exit 1; }
echo "== synthetic B/C pipeline test (fast; exercises mid/fine/C/fallback/final-C and pipeline_check)"
S2=$(date +%s.%N)
( cd $I/stageBC && bash build.sh brun.cpp "$I/stageBC/build/brun" && bash build.sh crun.cpp "$I/stageBC/build/crun" ) > /dev/null 2>&1 || { echo "SMOKE FAIL: stageBC rebuild"; exit 1; }
bash $I/stageBC/tests/pipeline_test.sh > $P/work_smoke_synth.log 2>&1 || { tail -5 $P/work_smoke_synth.log; echo "SMOKE FAIL: synthetic B/C pipeline test"; exit 1; }
tail -2 $P/work_smoke_synth.log
echo "[time] synthetic B/C: $(echo "$(date +%s.%N) - $S2" | bc) s"
echo "[time] smoke total: $(echo "$(date +%s.%N) - $S" | bc) s"
grep -q 'stage B/C INCOMPLETE' $L && echo "NOTE: B/C on the REAL Stage A records did not complete within the budget (see residual line above)"
echo "SMOKE OK (run_all rc=$rc)"
exit 0
