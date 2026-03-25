% basic scripts to start the resting state data analysis

restoredefaultpath; clear all
if ismac
    basepath = '/Users/kloosterman/Documents/GitHub/'; % local
    backend = 'none'; % local torque2    addpath(fullfile(basepath, 'MEG2afc'))
    addpath(genpath(fullfile(basepath, 'plotting-tools/')))
    addpath(genpath(fullfile(basepath, 'stats_tools/')))
    toolspath = '/Users/kloosterman/Documents/GitHub/';
    addpath(fullfile('/Users/kloosterman/Dropbox/tardis_code/MATLAB/tools/NoiseTools')) % robust detrend (MSE)
else
    basepath = '/mnt/beegfs/home/kloosterman/GitHub'; % on the cluster
%     addpath(fullfile(basepath, 'tools'))
    backend = 'slurm'; % local torque slurm
end
addpath(genpath(fullfile(basepath, 'entropyAge')))
addpath(fullfile(basepath, 'fieldtrip')) % cloned on 13 09 19
ft_defaults
addpath(fullfile(basepath, 'zapline-plus')) 
addpath(fullfile(basepath, 'qsub-tardis')) %inc JJ edit ft_artifact_zvalue

addpath('/Users/kloosterman/Documents/GitHub/fieldtrip_dev')
addpath(fullfile(toolspath, 'plscmd'))

set(groot,'defaultTextInterpreter','none')
set(groot,'defaultAxesTickLabelInterpreter','none')
set(groot,'defaultLegendInterpreter','none')
%% preprocessing EEG data
eA_preproc_setup() % done by Moritz

runPLSanalyses

%% plot YA vs OA taskPLS and time courses blink data
plot_taskPLS_YAvsOA 

%% plot behavior pls blink data

plot_PLSCblinkresults % behav PLS
plot_PLSCcorr_blink % corr bar plots
plotBS_taskPLSvsRawbehavior % task PLS BS vs behavior

