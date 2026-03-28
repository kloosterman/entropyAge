close all

folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

load young_mse_all.mat
load old_mse_all.mat
load behav.mat

% drop subjects if necessary


%% task PLS YA vs OA BLINK data
cfg = [];
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 100;                         % Number of permutation
cfg.num_boot = 100;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 1;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 0;                            % 0 is Pearson corr, 8 is Spearman
cfg.num_cond = 1;                           % Number of conditions
cfg.design = ones([1 39]);
cfg.num_subj_lst = [20 19];                 % Number of subjects per condition, array!
stat_mse_2group_task = ft_freqstatistics(cfg, young_mse_all, old_mse_all);

stat_mse_2group_task.posclusterslabelmat = stat_mse_2group_task.stat > 3 | stat_mse_2group_task.stat < -3;

% TODO save

%% Behavioral PLS analysis BLINK data
folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

out = []; % remove YA subject 2?
out_bool = true(20,1);
out_bool(out) = false;

loadstat = 0;
if loadstat == 1  
  load stat_mse_2group_blink.mat
  % load stat_mse_2group_blink_modulation.mat
else % run the behav pls here
  load young_mse_all.mat
  young_mse_all.powspctrm = young_mse_all.powspctrm(out_bool,:,:,:);
  load old_mse_all.mat

  % Split into groups
  behav_YA = behav(1:20,:);
  behav_OA = behav(21:end,:);

  out_bool_YA = true(size(behav_YA,1),1);
  out_bool_YA(out) = false;
  behav_YA = behav_YA(out_bool_YA,:);

  % (optional) OA removal placeholder
  % out_OA = []; % e.g. [3]
  % out_bool_OA = true(size(behav_OA,1),1);
  % out_bool_OA(out_OA) = false;
  % behav_OA = behav_OA(out_bool_OA,:);

  % Recombine
  behav = [behav_YA; behav_OA];

  cfg = [];
  cfg.frequency = [20 100];
  cfg.statistic = 'ft_statfun_pls';           % PLS statistics
  cfg.num_perm = 1000;                         % Number of permutation
  cfg.num_boot = 1000;
  cfg.method = 'analytic';                    % analytic method for statistics
  cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
  cfg.cormode = 0;                            % 0 is Pearson corr, 8 is Spearman
  cfg.num_cond = 1;                           % Number of conditions
  % cfg.design = behav;        % append behav OA to YA
  cfg.design = zscore(behav_domains);        % append behav OA to YA
  cfg.num_subj_lst = [size(young_mse_all.powspctrm,1) size(old_mse_all.powspctrm,1)];                 % Number of subjects per condition, array!
  stat_mse_2group_behav = ft_freqstatistics(cfg, young_mse_all, old_mse_all);

  stat_mse_2group_behav.posclusterslabelmat = stat_mse_2group_behav.stat > 3 | stat_mse_2group_behav.stat < -3;
end

% flip brain and behav scores YA for LV1 and 2
for ilv = 1:2
  stat_mse_2group_behav.brainscores{1}(:,ilv) = -stat_mse_2group_behav.brainscores{1}(:,ilv);
  stat_mse_2group_behav.behavscores{1}(:,ilv) = -stat_mse_2group_behav.behavscores{1}(:,ilv);
  % flip brain and behav scores OA
  stat_mse_2group_behav.brainscores{2}(:,ilv) = -stat_mse_2group_behav.brainscores{2}(:,ilv);
  stat_mse_2group_behav.behavscores{2}(:,ilv) = -stat_mse_2group_behav.behavscores{2}(:,ilv);
  % flip confidence intervals
  stat_mse_2group_behav.boot_res.ulcorr(:,ilv) = - stat_mse_2group_behav.boot_res.llcorr(:,ilv);
  stat_mse_2group_behav.results.boot_result.compare_u(:,ilv) = -stat_mse_2group_behav.results.boot_result.compare_u(:,ilv);
end

% flip BSRs too; LV1 by default in stat.stat
stat_mse_2group_behav.stat = -stat_mse_2group_behav.stat;

% tmp = stat_mse_2group_behav.boot_res.ulcorr;
% stat_mse_2group_behav.boot_res.ulcorr = - stat_mse_2group_behav.boot_res.llcorr;
% stat_mse_2group_behav.boot_res.llcorr = - tmp;

% save to disk
% save(fullfile(plotfolder, 'behavPLS_2group_blink'), 'stat_mse_2group_behav')

%% TODO no blink

%% TODO no blink modulation


