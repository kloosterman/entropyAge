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

%% compute cog domain and cog regime scores
compute_behav_domain_regime

%% plot mean ctr PLS and correlations with flex axis
plot_taskPLS_YAvsOA
plotBS_taskPLSvsRawbehavior % task PLS BS vs behavior
outfile = sprintf('taskPLS_corrtobehav');
exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), 'ContentType', 'vector', 'BackgroundColor','white')
exportgraphics(gcf, fullfile(plotfolder, [outfile '.png']), 'ContentType', 'vector', 'BackgroundColor','white')

%% plot behavior pls blink data
for LVsel = 1:3
  plot_PLSCblinkresults % behav PLS
  plot_PLSCcorr_blink % corr bar plots
  outfile = sprintf('LV%d_BehavPLS', LVsel);
  exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), 'ContentType', 'vector', 'BackgroundColor','white')
  exportgraphics(gcf, fullfile(plotfolder, [outfile '.png']), 'ContentType', 'vector', 'BackgroundColor','white')
end
