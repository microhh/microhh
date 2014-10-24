#!/bin/bash
#SBATCH --job-name=moser600
#SBATCH --output=mhh-%j.out
#SBATCH --error=mhh-%j.err
#SBATCH --partition=mpi-compute
#SBATCH --mail-type=end
#SBATCH --mail-user=abc@def.de
#SBATCH --nodes=32
#SBATCH -t 08:00:00 

# CASE SETTINGS
# ----------------------
case='moser600'
totaltime=5400      # total simulation time
maxtime=1800         # integration time per job
precision=1.00      
ncores=512 
# ----------------------

# 1. Check if .tmp file exist, if not copy case.ini
if [ ! -f $case.tmp ]; then
  echo "no $case.tmp file, copying $case.ini"
  cp $case.ini $case.tmp
fi

# 2. Get the latest time step of the restart files
timename=`ls time.* | tail -n1`
starttime="${timename##*.}"
starttime=`echo ${starttime} | sed 's/^0*//g'`  # remove 0000s
starttime=$(($starttime + 0))

# 3. Multiply with precision
digits=`echo "scale=0;l(1/${precision})/l(10)" | bc -l`
starttime=`echo "scale=${digits};${starttime}*${precision}" | bc`
runtime=`echo "scale=${digits};${starttime}+${maxtime}" | bc`

# 4. Limit the end time of the simulation at the total time
if [ `echo "${runtime}>${totaltime}" | bc` -eq "1" ]; then
  runtime=`echo "scale=${digits};${totaltime}" | bc`
fi

# 5. Create the ini file from the tmp file
echo "Inserting runtime ${runtime} into $case.ini"
sed "s/\(endtime\).*/\1=${runtime}/g" $case.tmp > $case.ini 
echo "Inserting starttime ${starttime} into $case.ini"
sed -i -e "s/\(starttime\).*/\1=${starttime}/g" $case.ini 

# 6. Init if t=0 + run the model
if [ ${starttime} -eq "0" ]; then
  mpirun --mca btl self,sm,openib -n $ncores ./microhh init $case
  mpirun --mca btl self,sm,openib -n $ncores ./microhh run $case
else
  mpirun --mca btl self,sm,openib -n $ncores ./microhh run $case
fi

# 7. Get return code model
runjobval=$?

# 8. Again get the latest time step of the restart files in case the simulation hit the wall clock limit
timename=`ls time.* | tail -n1`
starttime="${timename##*.}"
starttime=`echo ${starttime} | sed 's/^0*//g'`
starttime=$(($starttime + 0))

# If there is sim time left -> restart
timeleft=`echo "${starttime}<${totaltime}" | bc`
if [ "${runjobval}" -eq "0" -a "${timeleft}" -eq "1" ]; then
  sbatch run_restart.slurm
fi
