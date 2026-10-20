#!/bin/bash
# Synthetic end-to-end pipeline: Stage-A-like records (certified Roots at Fourier-family Hadamards, radius .004)
#   -> B mid -> B fine (batches) -> C per batch (tiny graph cap forces some UNRESOLVED) -> B fallback -> C final.
# Passes iff the run ends with "0 unresolved" or an explicit residual list, with every structural invariant intact.
set -u
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-8}
cd "$(dirname "$0")/.."
BIN=./build
W=build/pipe
rm -rf $W && mkdir -p $W
$BIN/mkrec fourier 0.3 0.5 0.004 $W/in0.bin > /dev/null || exit 1
$BIN/mkrec fourier -0.6 0.2 0.004 $W/in1.bin > /dev/null || exit 1
# the paper's parameters, except a small initial graph cap so that the fallback path is exercised
TAU_MID=.016 W_MID=.004 TAU_FINE=.012 W_FINE=.003 TAU_FB=.008 W_FB=.0025 CAP=${PIPE_CAP:-850} CAP_FINAL=5000000 CSHARDS=4 BIN=$BIN \
  bash ./run_pipeline.sh $W $W/in0.bin $W/in1.bin | tee $W/run.log
rc=${PIPESTATUS[0]}
python3 tests/pipeline_check.py $W .012 || exit 1
grep -q "^RESULT:" $W/run.log || { echo "pipeline produced no RESULT line"; exit 1; }
echo "initial-C unresolved batches: $(grep -c . <(grep 'batches with unresolved' $W/run.log))"
if grep -q "RESULT: 0 unresolved" $W/run.log; then echo "pipeline: 0 unresolved"; else echo "pipeline: residual list (rc=$rc)"; ls $W/RESIDUAL_* 2>/dev/null; fi
# the invariant "0 unresolved OR the residual list": rc must be 0 exactly when no RESIDUAL file exists
if [ $rc -eq 0 ]; then ls $W/RESIDUAL_* > /dev/null 2>&1 && { echo "rc=0 but residual files exist"; exit 1; }; else ls $W/RESIDUAL_* > /dev/null 2>&1 || { echo "rc!=0 but no residual file"; exit 1; }; fi
echo "pipeline_test: ok"
