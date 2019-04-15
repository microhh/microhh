#!/bin/bash

error_exit()
{
	echo "Failed" 1>&2
	exit 1
}

casename=${PWD##*/}
errorfile=$casename.err

echo "Case" $casename

$PYTHON_EXEC drycbllesprof.py >> $errorfile ||error_exit
rm -f *.000*
rm -f *.out
$MICROHH_EXEC init drycblles >> $errorfile ||error_exit
$MICROHH_EXEC run drycblles >> $errorfile ||error_exit
mv u.0003600 u.0003600ref
mv v.0003600 v.0003600ref
mv w.0003600 w.0003600ref
mv th.0003600 th.0003600ref
mv time.0003600 time.0003600ref
$MICROHH_EXEC run drycblles_restart >> $errorfile ||error_exit
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

