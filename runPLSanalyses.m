close all

folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

load young_mse_all.mat
load old_mse_all.mat
load behav.mat

% drop subjects if necessary


%% task PLS YA vs OA BLINK data
young_mse_all.dimord
size(young_mse_all.powspctrm)
size(old_mse_all.powspctrm)
% =========================================
% Optional subject exclusion
% ==========================================
ya_excl = [];          % e.g. [3 8]
oa_excl = [];         % e.g. [] for none

young_mse_use = young_mse_all;
old_mse_use   = old_mse_all;

if ~isempty(ya_excl)
    young_mse_use.powspctrm(ya_excl,:,:,:) = [];
end

if ~isempty(oa_excl)
    old_mse_use.powspctrm(oa_excl,:,:,:) = [];
end

nYA = size(young_mse_use.powspctrm,1);
nOA = size(old_mse_use.powspctrm,1);

cfg = [];
cfg.frequency    = [20 100];
cfg.statistic    = 'ft_statfun_pls';
cfg.num_perm     = 1000;
cfg.num_boot     = 1000;
cfg.method       = 'analytic';
cfg.pls_method   = 1;
cfg.cormode      = 0;
cfg.num_cond     = 1;
cfg.design       = ones(1, nYA + nOA);
cfg.num_subj_lst = [nYA nOA];

stat_mse_2group_task = ft_freqstatistics(cfg, young_mse_use, old_mse_use);

stat_mse_2group_task.posclusterslabelmat = ...
    stat_mse_2group_task.stat > 3 | stat_mse_2group_task.stat < -3;
stat_mse_2group_task_all = stat_mse_2group_task;

% now again without the outlier
% =========================================
% Optional subject exclusion
% ==========================================
ya_excl = [];          % e.g. [3 8]
oa_excl = 5;         % e.g. [] for none

young_mse_use = young_mse_all;
old_mse_use   = old_mse_all;

if ~isempty(ya_excl)
    young_mse_use.powspctrm(ya_excl,:,:,:) = [];
end

if ~isempty(oa_excl)
    old_mse_use.powspctrm(oa_excl,:,:,:) = [];
end

nYA = size(young_mse_use.powspctrm,1);
nOA = size(old_mse_use.powspctrm,1);

cfg = [];
cfg.frequency    = [20 100];
cfg.statistic    = 'ft_statfun_pls';
cfg.num_perm     = 1000;
cfg.num_boot     = 1000;
cfg.method       = 'analytic';
cfg.pls_method   = 1;
cfg.cormode      = 0;
cfg.num_cond     = 1;
cfg.design       = ones(1, nYA + nOA);
cfg.num_subj_lst = [nYA nOA];

stat_mse_2group_task = ft_freqstatistics(cfg, young_mse_use, old_mse_use);

stat_mse_2group_task.posclusterslabelmat = ...
    stat_mse_2group_task.stat > 3 | stat_mse_2group_task.stat < -3;
stat_mse_2group_task_noOAoutlier = stat_mse_2group_task;

% TODO save

%% Behavioral PLS analysis BLINK data
folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

% isOld = domain_z_tbl.AgeGroup == 'Older';
% cryst_OA = domain_z_tbl.Crystallized(isOld);
% idx_out = find(cryst_OA < -3);

out_YA = [];   % e.g. [2]
out_OA = 5;   % OA 5 has outlier low in Cr IQ

loadstat = 0;

if loadstat == 1
    load stat_mse_2group_blink.mat
    % load stat_mse_2group_blink_modulation.mat

else
    % Load brain data
    load young_mse_all.mat
    load old_mse_all.mat

    % Original group sizes
    nYA = size(young_mse_all.powspctrm, 1);
    nOA = size(old_mse_all.powspctrm, 1);

    % Inclusion masks
    keep_YA = true(nYA,1);
    keep_OA = true(nOA,1);
    keep_YA(out_YA) = false;
    keep_OA(out_OA) = false;

    % Apply exclusions to brain data
    young_mse_all.powspctrm = young_mse_all.powspctrm(keep_YA,:,:,:);
    old_mse_all.powspctrm   = old_mse_all.powspctrm(keep_OA,:,:,:);

    % Apply exclusions to behavior
    behav_YA = behav(1:nYA,:);
    behav_OA = behav(nYA+1:nYA+nOA,:);

    behav_YA = behav_YA(keep_YA,:);
    behav_OA = behav_OA(keep_OA,:);
    behav_clean = [behav_YA; behav_OA];

    % Apply exclusions to domain scores
    behav_domains_YA = behav_domains(1:nYA,:);
    behav_domains_OA = behav_domains(nYA+1:nYA+nOA,:);

    behav_domains_YA = behav_domains_YA(keep_YA,:);
    behav_domains_OA = behav_domains_OA(keep_OA,:);
    behav_domains_clean = [behav_domains_YA; behav_domains_OA];

    % Run behavioral PLS
    cfg = [];
    cfg.frequency    = [20 100];
    cfg.statistic    = 'ft_statfun_pls';
    cfg.num_perm     = 1000;
    cfg.num_boot     = 1000;
    cfg.method       = 'analytic';
    cfg.pls_method   = 3;   % 1 = taskPLS; 3 = behavPLS
    cfg.cormode      = 0;   % 0 = Pearson, 8 = Spearman
    cfg.num_cond     = 1;
    cfg.design       = behav_domains_clean;   % or behav_clean
    cfg.num_subj_lst = [sum(keep_YA), sum(keep_OA)];

    stat_mse_2group_behav = ft_freqstatistics(cfg, young_mse_all, old_mse_all);

    stat_mse_2group_behav.posclusterslabelmat = ...
        stat_mse_2group_behav.stat > 3 | stat_mse_2group_behav.stat < -3;
end

% Flip sign for interpretation
for ilv = [1 3] % 3 only flip when leaving outlier out
    % YA
    stat_mse_2group_behav.brainscores{1}(:,ilv) = -stat_mse_2group_behav.brainscores{1}(:,ilv);
    stat_mse_2group_behav.behavscores{1}(:,ilv) = -stat_mse_2group_behav.behavscores{1}(:,ilv);

    % OA
    stat_mse_2group_behav.brainscores{2}(:,ilv) = -stat_mse_2group_behav.brainscores{2}(:,ilv);
    stat_mse_2group_behav.behavscores{2}(:,ilv) = -stat_mse_2group_behav.behavscores{2}(:,ilv);

    % Flip confidence intervals correctly: swap + flip
    tmp_ul = stat_mse_2group_behav.boot_res.ulcorr(:,ilv);
    tmp_ll = stat_mse_2group_behav.boot_res.llcorr(:,ilv);

    stat_mse_2group_behav.boot_res.ulcorr(:,ilv) = -tmp_ll;
    stat_mse_2group_behav.boot_res.llcorr(:,ilv) = -tmp_ul;

    % Optional: flip compare intervals too, if lower bound exists
    if isfield(stat_mse_2group_behav.results,'boot_result')
        if isfield(stat_mse_2group_behav.results.boot_result,'compare_u') && ...
           isfield(stat_mse_2group_behav.results.boot_result,'compare_l')

            tmp_u = stat_mse_2group_behav.results.boot_result.compare_u(:,ilv);
            tmp_l = stat_mse_2group_behav.results.boot_result.compare_l(:,ilv);

            stat_mse_2group_behav.results.boot_result.compare_u(:,ilv) = -tmp_l;
            stat_mse_2group_behav.results.boot_result.compare_l(:,ilv) = -tmp_u;

        elseif isfield(stat_mse_2group_behav.results.boot_result,'compare_u')
            stat_mse_2group_behav.results.boot_result.compare_u(:,ilv) = ...
                -stat_mse_2group_behav.results.boot_result.compare_u(:,ilv);
        end
    end
end

% Flip BSR/stat map too
stat_mse_2group_behav.stat = -stat_mse_2group_behav.stat;

% Quick sanity check for CI order
for ilv = 1:2
    disp(all(stat_mse_2group_behav.boot_res.llcorr(:,ilv) <= ...
             stat_mse_2group_behav.boot_res.ulcorr(:,ilv)));
end

% save(fullfile(plotfolder, 'behavPLS_2group_blink.mat'), 'stat_mse_2group_behav')

%% TODO no blink

%% TODO no blink modulation


