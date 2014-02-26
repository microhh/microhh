#!/bin/bash
python drycbllesprof.py

rm -f *.000*
rm -f *.out
mpiexec -n 8 ./microhh init drycblles
mpiexec -n 8 ./microhh run drycblles
mv u.0003600 u.0003600ref
mv v.0003600 v.0003600ref
mv w.0003600 w.0003600ref
mv th.0003600 th.0003600ref
mv obuk.0003600 obuk.0003600ref
mpiexec -n 8 ./microhh run drycblles_restart
cmp u.0003600 u.0003600ref
diffu=$?
cmp v.0003600 v.0003600ref
diffv=$?
cmp w.0003600 w.0003600ref
diffw=$?
cmp th.0003600 th.0003600ref
diffs=$?
cmp obuk.0003600 obuk.0003600ref
diffobuk=$?
error=$(($diffu + $diffv + $diffw + $diffs + $diffobuk))
if [ $error = 0 ]; then
  echo "TEST PASSED!"
else
  echo "TEST FAILED!"
fi

