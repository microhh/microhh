#!/bin/bash
#SBATCH --job-name=rico
#SBATCH --output=mhh-%j.out
#SBATCH --error=mhh-%j.err
#SBATCH --partition=gpu
#SBATCH --gpus=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH -t 08:00:00

module purge
module load 2021
module load CMake/3.20.1-GCCcore-10.3.0
module load foss/2021a
module load netCDF/4.8.0-gompi-2021a
module load CUDA/11.3.1

./microhh init rico
./microhh run rico
