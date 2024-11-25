#!/bin/bash
#PBS -l walltime=47:00:00
#PBS -l select=1:mem=3gb


cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
source venv/bin/activate

python count_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-80-90-deuteron/"  3000 500 1000010020 
