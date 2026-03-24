close all

folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

load young_mse_all.mat
load old_mse_all.mat

load colormap_jetlightgray.mat

% cfg=[];
% cfg.layout = 'EEG1010.lay';
% cfg.interactive = 'yes';
% ft_multiplotTFR(cfg, young_mse_all)

%% YA vs OA stats task PLS
%% 2 group PLS
% older_mse_all = young_mse_all; % @Moritz: put the OA data here!
% older_mse_all.powspctrm = older_mse_all.powspctrm(1:17,:,:,:); % pretend lower N in OA

cfg = [];
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 100;                         % Number of permutation
cfg.num_boot = 100;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 1;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 0;                            % 0 is Pearson corr, 8 is Spearman
cfg.num_cond = 1;                           % Number of conditions
% cfg.interaction = 'yes'; % add group interaction to the model
% cfg.contrast = [-1 1];   % -1 for YA, 1 for OA
% cfg.design = [behav; behav(1:17,:)];        % append behav OA to YA
cfg.design = ones([1 39]);
cfg.num_subj_lst = [20 19];                 % Number of subjects per condition, array!
stat_mse_2group_task = ft_freqstatistics(cfg, young_mse_all, old_mse_all);

stat_mse_2group_task.posclusterslabelmat = stat_mse_2group_task.stat > 3 | stat_mse_2group_task.stat < -3;

%% plot 

close all

fig=figure;
tiledlayout(2,3)
fig.Position = [451   695   607   254];

cfg=[];
cfg.layout = 'EEG1010.lay';
cfg.clus2plot = 1;
cfg.clussign = 'pos';
cfg.integratetype = 'trapz'; % mean or trapz
cfg.subplotsize = [2 4]; % 2 rows, 2 topo/TFR couples
cfg.subplotind = [1 2];
cfg.parameter = 'stat';
cfg.colormap = cmap;
cfg.ylabel = 'Time scale (ms)';

ft_clusterplot3D(cfg, stat_mse_2group_task)

fig.Units = 'centimeters';
fig.Position(3:4) = [13 8];   % figure size: 12 × 8 cm

fig.PaperUnits = 'centimeters';
fig.PaperSize = [13 8];
fig.PaperPosition = [0 0 13 8];

exportgraphics(fig,fullfile(plotfolder, 'taskPLS_YAvsOA.pdf'),'ContentType','vector')

%% plot brain scores bar plot

Y_BStask = stat_mse_2group_task.brainscores{1}(:,1);
O_BStask = stat_mse_2group_task.brainscores{2}(:,1);

f=figure; 
tiledlayout(1,2)
nexttile
histogram(O_BStask,10)
hold on;
histogram(Y_BStask,10)
nexttile
bar([mean(O_BStask) mean(Y_BStask)])
