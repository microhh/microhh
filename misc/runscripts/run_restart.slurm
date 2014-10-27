#!/bin/bash
#SBATCH --job-name=abcjob
#SBATCH --output=mhh-%j.out
#SBATCH --error=mhh-%j.err
#SBATCH --partition=mpi-compute
#SBATCH --mail-type=end
#SBATCH --mail-user=abc@def.de
#SBATCH --nodes=32
#SBATCH -t 08:00:00

# CASE SETTINGS (EDIT HERE)
# ----------------------
case='thermal'
totaltime=2  # total simulation time
maxtime=0.5  # integration time per job
precision=0.1
nprocesses=512
runscriptname='run_restart.slurm'
# ----------------------

# 1. Copy .ini to .tmp
cp $case.ini $case.tmp

# 2. Get the latest time step of the restart files
timename=`ls time.* | tail -n1`
istarttime="${timename##*.}"
istarttime=`echo ${istarttime} | sed 's/^0*//g'`  # remove 0000s
istarttime=$(($istarttime + 0))

# 3. Multiply with precision
digits=`echo "scale=0;l(1/${precision})/l(10)" | bc -l`
starttime=`echo "scale=${digits};${istarttime}*${precision}" | bc`
runtime=`echo "scale=${digits};${starttime}+${maxtime}" | bc`

# 4. Limit the end time of the simulation at the total time
if [ `echo "${runtime}>${totaltime}" | bc` -eq "1" ]; then
  runtime=`echo "scale=${digits};${totaltime}" | bc`
fi

# 5. Create the ini file from the tmp file
echo "RUNSCRIPT: Inserting runtime ${runtime} into $case.ini"
sed "s/\(endtime\).*/\1=${runtime}/g" $case.tmp > $case.ini
echo "RUNSCRIPT: Inserting starttime ${starttime} into $case.ini"
sed -i -e "s/\(starttime\).*/\1=${starttime}/g" $case.ini

# 6. Init if t=0 + run the model
if [ ${istarttime} -eq "0" ]; then
  #mpirun --mca btl self,sm,openib -n $nprocesses ./microhh init $case
  #mpirun --mca btl self,sm,openib -n $nprocesses ./microhh run $case
  ./microhh init $case
  ./microhh run $case
else
  #mpirun --mca btl self,sm,openib -n $nprocesses ./microhh run $case
  ./microhh run $case
fi

# 7. Get return code model
runjobval=$?

# 8. Again get the latest time step of the restart files in case the simulation hit the wall clock limit
timename=`ls time.* | tail -n1`
istarttime="${timename##*.}"
istarttime=`echo ${istarttime} | sed 's/^0*//g'`
istarttime=$(($istarttime + 0))

# 9. Multiply with precision
digits=`echo "scale=0;l(1/${precision})/l(10)" | bc -l`
starttime=`echo "scale=${digits};${istarttime}*${precision}" | bc`

# If there is sim time left -> restart
timeleft=`echo "${starttime}<${totaltime}" | bc`
if [ "${runjobval}" -eq "0" -a "${timeleft}" -eq "1" ]; then
  echo "RUNSCRIPT: Dispatching job with starttime ${starttime}"
  #sbatch $runscriptname
  ./$runscriptname
else
  echo "RUNSCRIPT: Totaltime ${totaltime} has been reached."
fi

