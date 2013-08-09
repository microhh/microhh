#!/bin/bash
rm -f *.000*
rm -f *.out
mpiexec -n 4 ./microhh init drycblles
mpiexec -n 4 ./microhh run drycblles
mv u.0001800 u.0001800ref
mv v.0001800 v.0001800ref
mv w.0001800 w.0001800ref
mv s.0001800 s.0001800ref
mpiexec -n 4 ./microhh run drycblles_restart
cmp u.0001800 u.0001800ref
diffu=$?
cmp v.0001800 v.0001800ref
diffv=$?
cmp w.0001800 w.0001800ref
diffw=$?
cmp s.0001800 s.0001800ref
diffs=$?
error=$(($diffu + $diffv + $diffw + $diffs))
if [ $error = 0 ]; then
  echo "TEST PASSED!"
else
  echo "TEST FAILED!"
fi

