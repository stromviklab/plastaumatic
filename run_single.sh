#!/bin/bash
#SBATCH --account=rpp-mstrom
#SBATCH --time=06:00:00
#SBATCH --mem=128G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mail-user=wenyi.chen2@mail.mcgill.ca
#SBATCH --mail-type=ALL
#SBATCH --output=slurm_snakemake.out

module load StdEnv/2020
module load gcc/9.3.0
module load trimmomatic/0.39
module load quast/5.0.2
module load mummer/4.0.0beta2
module load samtools/1.13

time bash job_input.sh