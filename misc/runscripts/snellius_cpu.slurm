#!/bin/bash
#SBATCH --job-name=drycblles
#SBATCH --output=mhh-%j.out
#SBATCH --error=mhh-%j.err
#SBATCH --partition=thin
#SBATCH -n 128
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-core=1
#SBATCH -t 08:00:00
#SBATCH --mail-type=none
#SBATCH --mail-user=ceo@microhh.org

# Switch between GCC and Intel builds:
toolkit="gcc"     # gcc/intel

module purge
module load 2021

if [ "$toolkit" = "gcc" ]; then
    module load CMake/3.20.1-GCCcore-10.3.0
    module load foss/2021a
    module load netCDF/4.8.0-gompi-2021a
elif [ "$toolkit" = "intel" ]; then
    module load intel/2021a
    module load netCDF/4.8.0-iimpi-2021a
    module load FFTW/3.3.9-intel-2021a
fi

srun ./microhh init drycblles
srun ./microhh run drycblles
