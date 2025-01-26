#!/bin/bash
#PBS -l walltime=70:00:00
#PBS -l select=1:mem=27gb


cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
source venv/bin/activate

#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-00-005-deuteron/" "lhc2760-00-005" 3000 500 2 2212

python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-50-60-deuteron/" "lhc2760-50-60-python-1" 3000 500 2 2212