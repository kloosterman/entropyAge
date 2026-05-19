% Lightweight renderer for the manuscript figure assembled in eA_runanalysis.
% This skips preprocessing and recomputes only the task PLS variables needed
% by plot_taskPLS_slim.

restoredefaultpath
basepath = '/Users/kloosterman/Documents/GitHub/';
addpath(genpath(fullfile(basepath, 'plotting-tools/')))
addpath(genpath(fullfile(basepath, 'stats_tools/')))
addpath(genpath(fullfile(basepath, 'entropyAge')))
addpath(fullfile(basepath, 'fieldtrip'))
addpath(fullfile(basepath, 'fieldtrip_dev'))
addpath(fullfile(basepath, 'plscmd'))
ft_defaults

set(groot,'defaultTextInterpreter','none')
set(groot,'defaultAxesTickLabelInterpreter','none')
set(groot,'defaultLegendInterpreter','none')
cmap = cbrewer('div','RdBu',256);
cmap = flipud(cmap);

folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";
cd(folder)

load young_mse_all.mat
load old_mse_all.mat
load behav.mat

cachefile = fullfile(plotfolder, 'taskPLS_cache_all_noOAoutlier.mat');
if isfile(cachefile)
    load(cachefile, 'stat_mse_2group_task_all', 'stat_mse_2group_task_noOAoutlier')
else
    ya_excl = [];
    oa_excl = [];
    young_mse_use = young_mse_all;
    old_mse_use = old_mse_all;
    young_mse_use.powspctrm(ya_excl,:,:,:) = [];
    old_mse_use.powspctrm(oa_excl,:,:,:) = [];

    nYA = size(young_mse_use.powspctrm,1);
    nOA = size(old_mse_use.powspctrm,1);

    cfg = [];
    cfg.frequency = [20 100];
    cfg.statistic = 'ft_statfun_pls';
    cfg.num_perm = 1000;
    cfg.num_boot = 1000;
    cfg.method = 'analytic';
    cfg.pls_method = 1;
    cfg.cormode = 0;
    cfg.num_cond = 1;
    cfg.design = ones(1, nYA + nOA);
    cfg.num_subj_lst = [nYA nOA];

    stat_mse_2group_task = ft_freqstatistics(cfg, young_mse_use, old_mse_use);
    stat_mse_2group_task.posclusterslabelmat = ...
        stat_mse_2group_task.stat > 3 | stat_mse_2group_task.stat < -3;
    stat_mse_2group_task_all = stat_mse_2group_task;

    ya_excl = [];
    oa_excl = 5;
    young_mse_use = young_mse_all;
    old_mse_use = old_mse_all;
    young_mse_use.powspctrm(ya_excl,:,:,:) = [];
    old_mse_use.powspctrm(oa_excl,:,:,:) = [];

    nYA = size(young_mse_use.powspctrm,1);
    nOA = size(old_mse_use.powspctrm,1);
    cfg.design = ones(1, nYA + nOA);
    cfg.num_subj_lst = [nYA nOA];

    stat_mse_2group_task = ft_freqstatistics(cfg, young_mse_use, old_mse_use);
    stat_mse_2group_task.posclusterslabelmat = ...
        stat_mse_2group_task.stat > 3 | stat_mse_2group_task.stat < -3;
    stat_mse_2group_task_noOAoutlier = stat_mse_2group_task;

    save(cachefile, 'stat_mse_2group_task_all', 'stat_mse_2group_task_noOAoutlier')
end

compute_behav_domain_regime
close all

figure('Units','centimeters','Position',[5 5 8.5 13]);
tiledlayout(3,6,"TileSpacing","compact","Padding","tight")

stat_mse_2group_task = stat_mse_2group_task_all;

plot_behav_bars
plot_taskPLS_slim
plotWhat = 'regimes';
plotBS_taskPLSvsRawbehavior

outfile = 'behav_taskPLS_corrGradient_debug';
set(gcf, 'Color', 'w')
set(findall(gcf, 'type', 'axes'), 'Color', 'w')
exportFigure(gcf, fullfile(plotfolder, [outfile '.pdf']), 'Resolution', 600, 'FontName', 'Arial', 'FontSize', 8)
exportFigure(gcf, fullfile(plotfolder, [outfile '.png']), 'Resolution', 600, 'FontName', 'Arial', 'FontSize', 8)
