#!/bin/bash
#PBS -l walltime=47:00:00
#PBS -l select=1:mem=3gb


cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
source venv/bin/activate

python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-00-05-deuteron/" "lhc2760-00-05" 3000 500 2 321

