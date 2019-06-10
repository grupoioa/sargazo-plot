#!/bin/bash
#SBATCH -o graph.out
#SBATCH -e graph.err
#SBATCH -J test
#SBATCH -p workq2
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
module load herramientas/python/3.6

PATH=/home/mroldan/.conda/envs/carto/bin:$PATH

echo "Running program on $SLURM_CPUS_ON_NODE CPU cores"
#conda activate carto
srun python sar_plot.py zeta "ELEVACIÓN DEL MAR" ../sargazo_graph/lustre_test/Sargazo01_0001.nc "eleva" 5
#srun python multi.py
