#!/bin/bash
#PBS -l walltime=60:00:00
#PBS -l select=1:mem=15gb


cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
source venv/bin/activate

python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-05-10-deuteron/" "lhc2760-05-10-python" 3000 500 2 1000010020
#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-00-05-deuteron/" "lhc2760-00-05-python" 3000 500 2 1000010020
#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-10-20-deuteron/" "lhc2760-10-20-python" 3000 500 2 1000010020
#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-20-30-deuteron/" "lhc2760-20-30-python" 3000 500 2 1000010020
#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-30-40-deuteron/" "lhc2760-30-40-python" 3000 500 2 1000010020
#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-40-50-deuteron/" "lhc2760-40-50-python" 3000 500 2 1000010020
#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-50-60-deuteron/" "lhc2760-50-60-python" 3000 500 2 1000010020
#python vn_SP_pT_pid.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-60-70-deuteron/" "lhc2760-60-70-python" 3000 500 2 1000010020


