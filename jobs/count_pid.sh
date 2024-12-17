#!/bin/bash
#PBS -l walltime=47:00:00
#PBS -l select=1:mem=15gb


cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
source venv/bin/activate

python count_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-40-50-deuteron/" "count-lhc2760-40-50-deuteron"  3000 500 1000010020 
