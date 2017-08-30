#!/bin/sh

#  GaussProcess.sh
#SBATCH --output=GaussProcess.out
#SBATCH --job-name=GaussProcess
#SBATCH -t 12:00:00
#SBATCH --ntask=1
#SBATCH -c 5
#SBATCH --mem-per-cpu 10000
#SBATCH --exclude node017,node018
#SBATCH --mail-type=ALL
#SBATCH --mail-user=ehoseini@mit.edu

module add mit/matlab/2016b
matlab -nodisplay -signelCompThread -r "addpath(genpath('/home/ehoseini/MyCodes/')); \
addpath(genpath('/home/ehoseini/MyCodes/MatlabTools/'));\
data_folder='/om/user/ehoseini/MyData/[2017]BMMcourse/yumu_20160531_singleplane_83.5Hz_mika_processing/plane9';\
cd(data_folder);\
load('ForebrainGaus.mat');\
ForebrainGausEstim=ForebrainGausProcess.GausEstimate;\
ForebrainCorrMat=ForebrainGausProcess.CorrMatrixFreq;\
save('ForebrainGausEstim','ForebrainGausEstim','-v7.3');\
save('ForebrainCorrMat','ForebrainCorrMat','-v7.3');\
load('TactumGaus.mat');\
TactumGausEstim=TactumGaussProcess.GausEstimate;\
TactumCorrMat=TactumGaussProcess.CorrMatrixFreq;\
save('TactumGausEstim','TactumGausEstim','-v7.3');\
save('TactumCorrMat','TactumCorrMat','-v7.3');\
quit;"
