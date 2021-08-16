#!/usr/bin/bash 

#SBATCH --partition=bgmp        ### Partition 
#SBATCH --job-name=demux         ### Job Name 
#SBATCH --output=jobname%j    ### File in which to store job output  
#SBATCH --time=10:00:00       
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8 
#SBATCH --account=bgmp 

conda activate bgmp_py39

/usr/bin/time -v  ./demultiplex.py  
