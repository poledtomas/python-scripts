
#PBS -l walltime=60:00:00
#PBS -l select=1:mem=3gb

module add root
cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
root -q -b 'neco.C("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-40-50-deuteron/",3000,500,2,1000010020,"lhc2760-40-50-root2")'