#!/bin/bash -l
# The -l above is required to get the full environment with modules

#SBATCH -A snic2022-3-25
#SBATCH  -p shared

# The name of the script is myjob
#SBATCH -J neko_post_grad

# Only 1 hour wall-clock time will be given to this job
#SBATCH -t 00:05:00

# Number of nodes
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
# Number of MPI processes.
#SBATCH -n 128
#SBATCH --exclusive

#SBATCH -e error_file.e
#SBATCH -o output_file.o
#SBATCH --mail-type=ALL
#SBATCH --mail-user=sdai@kth.se

# change modules
module swap PrgEnv-cray PrgEnv-gnu
module load PDC

casename='neko_post_grad'
numproc=128
echo $casename > SESSION.NAME
echo `pwd`'/' >> SESSION.NAME
# Run the executable
srun -n $numproc ./postprocess_fluid_stats box.nmsh mean_field1101_avg_z_x0.fld stats1101_avg_z_x0.fld mean > $casename.log 2>&1
