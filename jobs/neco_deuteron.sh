
#PBS -l walltime=60:00:00
#PBS -l select=1:mem=10gb

module add root
cd /storage/brno12-cerit/home/poledto1/hydro/hybrid/root_scripts/python-scripts
root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-00-05-deuteron/",3000,500,2,1000010020,"lhc2760-00-05")'
#root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-05-10-deuteron/",3000,500,2,1000010020,"lhc2760-05-10")'
#root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-10-20-deuteron/",3000,500,2,1000010020,"lhc2760-10-20")'
#root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-20-30-deuteron/",3000,500,2,1000010020,"lhc2760-20-30")'
#root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-30-40-deuteron/",3000,500,2,1000010020,"lhc2760-30-40")'
#root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-40-50-deuteron/",3000,500,2,1000010020,"lhc2760-40-50")'
#root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-50-60-deuteron/",3000,500,2,1000010020,"lhc2760-50-60")'
#root -q -b 'vn_SP_pT_pid.c("/storage/brno12-cerit/home/poledto1/hydro/hybrid/sampler.out/lhc2760-60-70-deuteron/",3000,500,2,1000010020,"lhc2760-60-70")'