#!/bin/bash

error_exit()
{
	echo "Failed" 1>&2
	exit 1
}

casename=${PWD##*/}
errorfile=$casename.err

echo "Case" $casename

cd conservation100_3rd
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..

cd conservation200_3rd
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..

cd conservation400_3rd
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..

cd conservation800_3rd
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..

cd conservation100_4th
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..

cd conservation200_4th
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..

cd conservation400_4th
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..

cd conservation800_4th
$PYTHON_EXEC conservationprof.py  >> $errorfile ||error_exit
rm -f *.00* conservation.out
$MICROHH_EXEC init conservation >> $errorfile ||error_exit
$MICROHH_EXEC run conservation >> $errorfile ||error_exit
cd ..
