#!/bin/bash
#PBS -l walltime=47:00:00
#PBS -l select=1:mem=3gb

cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts/
source venv/bin/activate

python spectrum.py "/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-30-40-deuteron/" "lhc2760-30-40-3000-deuteron-test" 3000 500