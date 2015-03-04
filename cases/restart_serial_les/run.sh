#!/bin/bash
python drycbllesprof.py
rm -f *.000*
rm -f *.out
./microhh init drycblles
./microhh run drycblles
mv u.0003600 u.0003600ref
mv v.0003600 v.0003600ref
mv w.0003600 w.0003600ref
mv th.0003600 th.0003600ref
mv time.0003600 time.0003600ref
./microhh run drycblles_restart
cmp u.0003600 u.0003600ref
diffu=$?
cmp v.0003600 v.0003600ref
diffv=$?
cmp w.0003600 w.0003600ref
diffw=$?
cmp th.0003600 th.0003600ref
diffs=$?
error=$(($diffu + $diffv + $diffw + $diffs))
if [ $error = 0 ]; then
  echo "TEST PASSED!"
else
  echo "TEST FAILED!"
fi

