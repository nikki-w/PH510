#!/bin/bash

#======================================================
#
# Job script for running a serial job on a single core 
#
#======================================================

#======================================================
# Propogate environment variables to the compute node
#SBATCH --export=ALL
#
# Run in the standard partition (queue)
#SBATCH --partition=teaching-gpu
#
# Specify project account
#SBATCH --account=teaching
#
# No. of tasks required (ntasks=1 for a single-core job)
#SBATCH --ntasks=16 --distribution=block:block
#
# Specify (hard) runtime (HH:MM:SS)
#SBATCH --time=00:30:00  
#
# Job name
#SBATCH --job-name=assignment1
#
# Output file
#SBATCH --output=slurm-%j.out
#======================================================

module purge
#Example module load command. 
#Load any modules appropriate for your program's requirements
module load openmpi/gcc-8.5.0/4.1.1
module add miniconda/3.12.8


#======================================================
# Prologue script to record job details
# Do not change the line below
#======================================================
/opt/software/scripts/job_prologue.sh  
#------------------------------------------------------

# Modify the line below to run your program
mpiexec -n 16 python3 ./assignment1.py

#======================================================
# Epilogue script to record job endtime and runtime
# Do not change the line below
#======================================================
/opt/software/scripts/job_epilogue.sh 
------------------------------------------------------
