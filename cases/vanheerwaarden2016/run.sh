#!/bin/bash

error_exit()
{
	echo "Failed" 1>&2
	exit 1
}

for dir in */
do
  cd $dir
  casename=${PWD##*/}
  errorfile=$casename.err

  echo "Case" $casename
  echo "Python preprocessing..."
  $PYTHON_EXEC ${casename}prof.py >> $errorfile ||error_exit
  echo "Model initialization..."
  ${MICROHH_EXEC} init thermal >> $errorfile ||error_exit
  echo "Model execution..."
  ${MICROHH_EXEC} run thermal >> $errorfile ||error_exit
  cd ..
done

echo "Succes!"

