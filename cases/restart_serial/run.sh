#!/bin/bash

error_exit()
{
	echo "Failed" 1>&2
	exit 1
}

casename=${PWD##*/}
errorfile=$casename.err

echo "Case" $casename

$PYTHON_EXEC drycbl_flowprof.py >> $errorfile ||error_exit

rm -f *.000*
rm -f *.out
$MICROHH_EXEC init drycbl_flow >> $errorfile ||error_exit
$MICROHH_EXEC run drycbl_flow >> $errorfile ||error_exit
mv u.0000002 u.0000002ref
mv v.0000002 v.0000002ref
mv w.0000002 w.0000002ref
mv b.0000002 b.0000002ref
mv time.0000002 time.0000002ref
$MICROHH_EXEC run drycbl_flow_restart >> $errorfile ||error_exit
cmp u.0000002 u.0000002ref
diffu=$?
cmp v.0000002 v.0000002ref
diffv=$?
cmp w.0000002 w.0000002ref
diffw=$?
cmp b.0000002 b.0000002ref
diffs=$?
error=$(($diffu + $diffv + $diffw + $diffs))
if [ $error = 0 ]; then
  echo "TEST PASSED!"
else
  echo "TEST FAILED!"
fi

