#!/bin/bash
#PBS -l walltime=60:00:00
#PBS -l select=1:mem=8gb

module add python
cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
source venv/bin/activate
python vn_EP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-30-40-deuteron" "lhc2760-30-40-3000-deuteron" 3000 500 2 1000010020 