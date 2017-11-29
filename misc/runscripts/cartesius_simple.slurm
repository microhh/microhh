#!/bin/bash
#SBATCH --job-name=casename
#SBATCH --output=mhh-%j.out
#SBATCH --error=mhh-%j.err
#SBATCH --partition=short
#SBATCH --mail-type=end
#SBATCH --mail-user=abc@def.de
#SBATCH -N 1
#SBATCH --ntasks-per-node=24
#SBATCH -t 01:00:00

# For now needed at Cartesius
export HDF5_DISABLE_VERSION_CHECK=1

# Get the case name from the .ini file
inifile=`ls *.ini | tail -n1`
case=$(echo $inifile | tr "." "\n")

# Run the model in init mode (only for t=0)
srun ./microhh init $case

# Run the model
srun ./microhh run $case
