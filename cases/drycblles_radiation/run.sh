#!/bin/bash

error_exit()
{
	echo "Failed" 1>&2
	exit 1
}

casename=drycblles
errorfile=$casename.err

echo "Case" $casename
echo "Python preprocessing..."
$PYTHON_EXEC ${casename}prof.py >> $errorfile ||error_exit
echo "Model initialization..."
${MICROHH_EXEC} init $casename >> $errorfile ||error_exit
echo "Model execution..."
${MICROHH_EXEC} run $casename >> $errorfile ||error_exit
echo "Succes!"

