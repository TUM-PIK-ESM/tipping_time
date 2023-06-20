#!/bin/bash 
#SBATCH --qos=priority 
#SBATCH --job-name=job2
#SBATCH --account=tipes
#SBATCH --mem=10000
#SBATCH --output=job-%j.log
#SBATCH --error=job-%j.err
#SBATCH --workdir=/home/mayayami/TT_fingerprints/


# sh cdo_script.sh
python TT_get_fingerprints.py