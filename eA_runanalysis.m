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

%% plotting scripts below
% The big figure with behavior, task PLS brain scores, gradient and example
% scatters
compute_behav_domain_regime % plot behav bar plot
plot_behav_PLS_grad_scat % plot task PLS bar + topo

%% better: behavior bars, task pls, gradient plot: stable to flexible
close all
figure('Units','centimeters','Position',[5 5 8.5 13]);
tiledlayout(3,6,"TileSpacing","compact","Padding","tight")

% stat_mse_2group_task = stat_mse_2group_task_noOAoutlier;
stat_mse_2group_task = stat_mse_2group_task_all; % include all

plot_behav_bars % panel A
plot_taskPLS_slim
% plot_flexgradient_BSvsCog
plotWhat = 'regimes';
plotBS_taskPLSvsRawbehavior

outfile = sprintf('behav_taskPLS_corrGradient');
set(gcf, 'Color', 'w')                    % figure background
set(findall(gcf, 'tygcpe', 'axes'), 'Color', 'w')  % all axes
exportFigure(gcf, fullfile(plotfolder, [outfile '.pdf']), 'Resolution', 600, 'FontName', 'Arial', 'FontSize', 8)
exportFigure(gcf, fullfile(plotfolder, [outfile '.png']), 'Resolution', 600, 'FontName', 'Arial', 'FontSize', 8)

%% plot taskBS vs behavior per domain
% close all
plotWhat = 'domains';
% plotWhat = 'regimes';
figure
plotBS_taskPLSvsRawbehavior


%% plot YA vs OA taskPLS and time courses blink data
plot_taskPLS_YAvsOA 

%% compute cog domain and cog regime scores
compute_behav_domain_regime


%% plot mean ctr PLS and correlations with flex axis
plot_taskPLS_YAvsOA
plotBS_taskPLSvsRawbehavior % task PLS BS vs behavior
outfile = sprintf('taskPLS_corrtobehav');
set(gcf, 'Color', 'w')                    % figure background
set(findall(gcf, 'tygcpe', 'axes'), 'Color', 'w')  % all axes
exportFigure(gcf, fullfile(plotfolder, [outfile '.pdf']), 'Resolution', 600, 'FontName', 'Arial', 'FontSize', 8)
exportFigure(gcf, fullfile(plotfolder, [outfile '.png']), 'Resolution', 600, 'FontName', 'Arial', 'FontSize', 8)
% exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), 'ContentType', 'vector', 'BackgroundColor','white')
% exportgraphics(gcf, fullfile(plotfolder, [outfile '.png']), 'ContentType', 'vector', 'BackgroundColor','white')

%% plot behavior pls blink data
for LVsel = 1:3 %1:4
  close all
  plot_PLSCblinkresults % behav PLS
  plot_PLSCcorr_blink % corr bar plots
end

