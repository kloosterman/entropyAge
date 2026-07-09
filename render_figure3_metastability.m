%% Figure 3: blink-locked entropy dynamics and metastability summary
%
% Builds a three-panel figure:
%   A) anterior/posterior entropy time courses with group SEM
%   B) LV1-LV3 pre/quench/recovery state switch
%   C) compact metastability-relevant metrics

close all

codefolder = "/Users/kloosterman/Documents/GitHub/entropyAge";
datafolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
noblinkfolder = "/Users/kloosterman/projectdata/EntropyAging/usedinpaper";
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

if exist(codefolder, 'dir')
    addpath(codefolder);
end
if ~exist(plotfolder, 'dir')
    mkdir(plotfolder);
end

load(fullfile(datafolder, "young_mse_all.mat"), "young_mse_all");
load(fullfile(datafolder, "old_mse_all.mat"), "old_mse_all");
S = load(fullfile(noblinkfolder, "young_mse_all_noblink.mat"), "young_mse_all");
young_mse_all_noblink = S.young_mse_all;
S = load(fullfile(noblinkfolder, "old_mse_all_noblink.mat"), "old_mse_all");
old_mse_all_noblink = S.old_mse_all;
stat_mse_2group_behav = load_behav_pls_stat(datafolder, plotfolder);

col_YA = [0.85 0.25 0.25];
col_OA = [0.20 0.45 0.85];
col_YA_light = lighten_color(col_YA, 0.58);
col_OA_light = lighten_color(col_OA, 0.58);

anterior_labels = {'Fp1','Fp2','AF7','AF3','AFz','AF4','AF8', ...
    'F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
    'FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8', ...
    'T7','C5','C3','C1','Cz','C2','C4','C6','T8'};
posterior_labels = {'TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8', ...
    'P7','P5','P3','P1','Pz','P2','P4','P6','P8', ...
    'PO7','PO3','POz','PO4','PO8','O1','Oz','O2'};

time = young_mse_all.time(:)';

% Panel A follows the existing exploratory code: anterior slow scales and
% posterior fast scale.
tc_ant_YA = extract_roi_timecourse(young_mse_all, anterior_labels, [40 100]);
tc_ant_OA = extract_roi_timecourse(old_mse_all, anterior_labels, [40 100]);
tc_post_YA = extract_roi_timecourse(young_mse_all, posterior_labels, 20);
tc_post_OA = extract_roi_timecourse(old_mse_all, posterior_labels, 20);

ref_ant_YA = extract_roi_reference(young_mse_all_noblink, anterior_labels, [40 100]);
ref_ant_OA = extract_roi_reference(old_mse_all_noblink, anterior_labels, [40 100]);
ref_post_YA = extract_roi_reference(young_mse_all_noblink, posterior_labels, 20);
ref_post_OA = extract_roi_reference(old_mse_all_noblink, posterior_labels, 20);

tc_ant_YA_psc = percent_signal_change(tc_ant_YA, ref_ant_YA);
tc_ant_OA_psc = percent_signal_change(tc_ant_OA, ref_ant_OA);
tc_post_YA_psc = percent_signal_change(tc_post_YA, ref_post_YA);
tc_post_OA_psc = percent_signal_change(tc_post_OA, ref_post_OA);

% Panel B: project all-channel x timescale snapshots onto fixed LV1/LV3
% axes made by time-collapsing the behavior-PLS saliences.
state_space_labels = stat_mse_2group_behav.label;
[traj_YA, traj_OA, state_info] = compute_state_space( ...
    young_mse_all, old_mse_all, stat_mse_2group_behav, state_space_labels, [20 100]);

metric_window = [-1.00 1.00];

rise_base = [-1.00 -0.55];
rise_pre = [-0.50 -0.05];
quench_pre = [-0.25 0.00];
quench_post = [0.10 0.40];
recovery_win = [0.50 1.00];
switch_windows = [quench_pre; quench_post; recovery_win];
reset_win = quench_post;

rise_YA = window_mean(tc_ant_YA, time, rise_pre) - window_mean(tc_ant_YA, time, rise_base);
rise_OA = window_mean(tc_ant_OA, time, rise_pre) - window_mean(tc_ant_OA, time, rise_base);

quench_YA = window_mean(tc_ant_YA_psc, time, quench_pre) - window_min(tc_ant_YA_psc, time, quench_post);
quench_OA = window_mean(tc_ant_OA_psc, time, quench_pre) - window_min(tc_ant_OA_psc, time, quench_post);

lv3_reset_YA = compute_lv3_reset_amplitude(traj_YA, time, rise_base, reset_win);
lv3_reset_OA = compute_lv3_reset_amplitude(traj_OA, time, rise_base, reset_win);
speed_peak_YA = compute_peak_state_speed(traj_YA, time, reset_win);
speed_peak_OA = compute_peak_state_speed(traj_OA, time, reset_win);
outside_time_YA = compute_time_outside_baseline_basin(traj_YA, time, rise_base, reset_win);
outside_time_OA = compute_time_outside_baseline_basin(traj_OA, time, rise_base, reset_win);
[centroid_stats, centroid_data] = compute_centroid_permutation(traj_YA, traj_OA, time, metric_window, 20000);
[occupation_stats, occupation_data] = compute_occupation_area_permutation( ...
    traj_YA, traj_OA, time, metric_window, 20000);
state_entropy_stats = compute_state_space_entropy_permutation( ...
    traj_YA, traj_OA, time, metric_window, 4:9, 20000);
transition_entropy_stats = compute_state_transition_entropy_permutation( ...
    traj_YA, traj_OA, time, metric_window, 4:9, 20000);
[escape_return_stats, escape_return_data] = compute_escape_return_permutation( ...
    traj_YA, traj_OA, time, rise_base, quench_post, [0.10 1.00], 20000);
[switch_states_YA, switch_metrics_YA] = compute_state_switch(traj_YA, time, switch_windows);
[switch_states_OA, switch_metrics_OA] = compute_state_switch(traj_OA, time, switch_windows);

stats_table = make_metric_table( ...
    ["LV3 reset amplitude"; "Peak state-space speed"; "Time outside baseline basin"; ...
     "Blink-free PSC min quench"; "Pre-blink entropy rise"], ...
    {lv3_reset_YA, speed_peak_YA, outside_time_YA, quench_YA, rise_YA}, ...
    {lv3_reset_OA, speed_peak_OA, outside_time_OA, quench_OA, rise_OA});
switch_stats_table = make_metric_table( ...
    ["LV1 state switch"; "LV3 state switch"; "2-D switch magnitude"; "Recovery distance"; ...
     "Pre LV1 state"; "Quench LV1 state"; "Recovery LV1 state"; ...
     "Pre LV3 state"; "Quench LV3 state"; "Recovery LV3 state"], ...
    {switch_metrics_YA.lv1_switch, switch_metrics_YA.lv3_switch, ...
     switch_metrics_YA.switch_magnitude, switch_metrics_YA.recovery_distance, ...
     switch_states_YA(:,1,1), switch_states_YA(:,2,1), switch_states_YA(:,3,1), ...
     switch_states_YA(:,1,2), switch_states_YA(:,2,2), switch_states_YA(:,3,2)}, ...
    {switch_metrics_OA.lv1_switch, switch_metrics_OA.lv3_switch, ...
     switch_metrics_OA.switch_magnitude, switch_metrics_OA.recovery_distance, ...
     switch_states_OA(:,1,1), switch_states_OA(:,2,1), switch_states_OA(:,3,1), ...
     switch_states_OA(:,1,2), switch_states_OA(:,2,2), switch_states_OA(:,3,2)});
writetable(stats_table, fullfile(plotfolder, "Figure3_metastability_metrics.csv"));
writetable(switch_stats_table, fullfile(plotfolder, "Figure3_state_switch_metrics.csv"));
writetable(centroid_stats, fullfile(plotfolder, "Figure3_state_space_centroid_test.csv"));
writetable(occupation_stats, fullfile(plotfolder, "Figure3_state_space_occupation_test.csv"));
writetable(state_entropy_stats, fullfile(plotfolder, "Figure3_state_space_entropy_test.csv"));
writetable(transition_entropy_stats, fullfile(plotfolder, "Figure3_state_transition_entropy_test.csv"));
writetable(escape_return_stats, fullfile(plotfolder, "Figure3_escape_return_test.csv"));

%% Figure layout

f = figure('Color', 'w');
f.Units = 'centimeters';
f.Position = [2 2 18.5 14.2];

axA1 = axes(f, 'Position', [0.075 0.635 0.385 0.295]);
plot_entropy_timecourse(axA1, time, tc_ant_YA_psc, tc_ant_OA_psc, col_YA, col_OA, ...
    'A  Anterior entropy', 'Change from blink-free (%)', rise_base, rise_pre, quench_pre, quench_post);

axA2 = axes(f, 'Position', [0.555 0.635 0.385 0.295]);
plot_entropy_timecourse(axA2, time, tc_post_YA_psc, tc_post_OA_psc, col_YA, col_OA, ...
    'Posterior entropy', 'Change from blink-free (%)', rise_base, rise_pre, quench_pre, quench_post);

yl = mean_sem_ylim({tc_ant_YA_psc, tc_ant_OA_psc, tc_post_YA_psc, tc_post_OA_psc}, 0.08);
axA1.YLim = yl;
axA2.YLim = yl;

axB = axes(f, 'Position', [0.075 0.110 0.500 0.405]);
plot_state_switch(axB, switch_states_YA, switch_states_OA, col_YA, col_OA);
title(axB, 'B  Blink-evoked LV state switch', 'FontWeight', 'bold');

axC1 = axes(f, 'Position', [0.680 0.415 0.255 0.105]);
plot_metric_panel(axC1, quench_YA, quench_OA, col_YA, col_OA, ...
    'C  Anterior quench', 'PSC drop (%)', stats_table.P_Group(4));

axC2 = axes(f, 'Position', [0.680 0.260 0.255 0.105]);
plot_metric_panel(axC2, switch_states_YA(:,2,1), switch_states_OA(:,2,1), col_YA, col_OA, ...
    'Quench LV1 state', 'LV1 expression', switch_stats_table.P_Group(6));

axC3 = axes(f, 'Position', [0.680 0.105 0.255 0.105]);
plot_metric_panel(axC3, switch_states_YA(:,3,1), switch_states_OA(:,3,1), col_YA, col_OA, ...
    'Recovery LV1 state', 'LV1 expression', switch_stats_table.P_Group(7));

annotation(f, 'textbox', [0.680 0.540 0.255 0.035], ...
    'String', sprintf('LV axes: %s', state_info.axis_label), ...
    'EdgeColor', 'none', 'FontSize', 6, 'Color', [0.25 0.25 0.25], ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

outfile = "Figure3_metastability";
if exist('exportFigure', 'file')
    exportFigure(f, fullfile(plotfolder, outfile + ".pdf"), 'Resolution', 600, ...
        'FontName', 'Arial', 'FontSize', 8);
    exportFigure(f, fullfile(plotfolder, outfile + ".png"), 'Resolution', 600, ...
        'FontName', 'Arial', 'FontSize', 8);
else
    exportgraphics(f, fullfile(plotfolder, outfile + ".pdf"), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f, fullfile(plotfolder, outfile + ".png"), ...
        'Resolution', 600, 'BackgroundColor', 'white');
end

fprintf('Saved %s.pdf/png and Figure3_metastability_metrics.csv to %s\n', ...
    outfile, char(plotfolder));

%% Local helpers

function stat = load_behav_pls_stat(datafolder, plotfolder)
candidate_files = [
    fullfile(plotfolder, "behavPLS_2group_blink.mat")
    fullfile(datafolder, "behavPLS_2group_blink.mat")
    fullfile(datafolder, "stats_Moritz", "stat_mse_2group_blink.mat")
    fullfile(datafolder, "stats_Moritz", "stat_mse_2group_blink_modulation.mat")
    ];

loaded = false;
for ifile = 1:numel(candidate_files)
    if exist(candidate_files(ifile), 'file')
        S = load(candidate_files(ifile));
        if isfield(S, 'stat_mse_2group_behav')
            stat = S.stat_mse_2group_behav;
            loaded = true;
            break
        elseif isfield(S, 'stat_mse_2group')
            stat = S.stat_mse_2group;
            loaded = true;
            break
        end
    end
end

if ~loaded
    error('Could not find a behavioral PLS stat struct.');
end
end

function tc = extract_roi_timecourse(data, roi_labels, freqrange)
chan_idx = ismember(data.label, roi_labels);
if ~any(chan_idx)
    error('No requested ROI labels were found.');
end

if isscalar(freqrange)
    freq_idx = abs(data.freq - freqrange) == min(abs(data.freq - freqrange));
else
    freq_idx = data.freq >= freqrange(1) & data.freq <= freqrange(2);
end
if ~any(freq_idx)
    error('No frequency/time-scale bins found for requested range.');
end

x = data.powspctrm(:, chan_idx, freq_idx, :);
tc = squeeze(mean(mean(x, 3, 'omitnan'), 2, 'omitnan'));
if isvector(tc)
    tc = tc(:);
end
end

function ref = extract_roi_reference(data, roi_labels, freqrange)
tc = extract_roi_timecourse(data, roi_labels, freqrange);
ref = mean(tc, 2, 'omitnan');
end

function tc_psc = percent_signal_change(tc, ref)
if size(tc, 1) ~= numel(ref)
    error('Blink-locked time courses and no-blink reference have different subject counts.');
end
ref = ref(:);
tc_psc = 100 .* (tc - ref) ./ ref;
end

function [traj_YA, traj_OA, info] = compute_state_space(dataYA, dataOA, stat, roi_labels, scale_ranges)
statdims = size(stat.stat);
if numel(statdims) ~= 3
    error('Expected stat.stat to be chan x freq x time.');
end
if size(scale_ranges, 1) == 1
    scale_ranges = repmat(scale_ranges, 2, 1);
end
if size(scale_ranges, 1) ~= 2 || size(scale_ranges, 2) ~= 2
    error('Expected scale_ranges to be 1x2 or 2x2.');
end

[tfY, data_chan_idx_YA] = ismember(stat.label, dataYA.label);
[tfO, data_chan_idx_OA] = ismember(stat.label, dataOA.label);
if ~all(tfY) || ~all(tfO)
    error('Could not align all PLS labels to the entropy data labels.');
end

data_freq_idx_YA = match_numeric_axis(stat.freq, dataYA.freq);
data_freq_idx_OA = match_numeric_axis(stat.freq, dataOA.freq);
roi_idx = ismember(stat.label, roi_labels);
time_axis_idx = stat.time >= -1.00 & stat.time <= 0.75;

if ~any(roi_idx) || ~any(time_axis_idx)
    error('Empty ROI, frequency, or time selection for state-space projection.');
end

nYA = size(dataYA.powspctrm, 1);
nOA = size(dataOA.powspctrm, 1);
nTime = size(dataYA.powspctrm, 4);
traj_YA = nan(nYA, nTime, 2);
traj_OA = nan(nOA, nTime, 2);

for ilv = 1:2
    lv = [1 3];
    freq_idx = stat.freq >= scale_ranges(ilv, 1) & stat.freq <= scale_ranges(ilv, 2);
    if ~any(freq_idx)
        error('Empty time-scale selection %.0f-%.0f ms for LV%d.', ...
            scale_ranges(ilv, 1), scale_ranges(ilv, 2), lv(ilv));
    end

    YA = dataYA.powspctrm(:, data_chan_idx_YA(roi_idx), data_freq_idx_YA(freq_idx), :);
    OA = dataOA.powspctrm(:, data_chan_idx_OA(roi_idx), data_freq_idx_OA(freq_idx), :);
    nFeat = nnz(roi_idx) * nnz(freq_idx);

    all_data = cat(1, YA, OA);
    Xall = reshape(permute(all_data, [1 4 2 3]), [], nFeat);
    mu = mean(Xall, 1, 'omitnan');
    sd = std(Xall, 0, 1, 'omitnan');
    sd(sd == 0 | isnan(sd)) = 1;

    XYA = (reshape(permute(YA, [1 4 2 3]), [], nFeat) - mu) ./ sd;
    XOA = (reshape(permute(OA, [1 4 2 3]), [], nFeat) - mu) ./ sd;

    u = reshape(stat.results.u(:, lv(ilv)), statdims);
    w = squeeze(mean(u(roi_idx, freq_idx, time_axis_idx), 3, 'omitnan'));
    w = w(:);
    if norm(w) == 0 || all(isnan(w))
        error('LV%d state-space axis has zero norm.', lv(ilv));
    end
    w = w ./ norm(w);

    proj_YA = reshape(XYA * w, nYA, nTime);
    proj_OA = reshape(XOA * w, nOA, nTime);

    all_proj = [proj_YA(:); proj_OA(:)];
    proj_mu = mean(all_proj, 'omitnan');
    proj_sd = std(all_proj, 0, 'omitnan');
    if proj_sd == 0 || isnan(proj_sd)
        proj_sd = 1;
    end

    traj_YA(:,:,ilv) = (proj_YA - proj_mu) ./ proj_sd;
    traj_OA(:,:,ilv) = (proj_OA - proj_mu) ./ proj_sd;
end

% Orient LV3 so the older peri-blink excursion is positive in the plot.
base_idx = stat.time >= -1.00 & stat.time <= -0.55;
reset_idx = stat.time >= -0.20 & stat.time <= 0.40;
oa_reset = mean(traj_OA(:, reset_idx, 2), 'all', 'omitnan') - ...
    mean(traj_OA(:, base_idx, 2), 'all', 'omitnan');
if oa_reset < 0
    traj_YA(:,:,2) = -traj_YA(:,:,2);
    traj_OA(:,:,2) = -traj_OA(:,:,2);
end

info.axis_label = sprintf('LV1 %.0f-%.0f ms; LV3 %.0f-%.0f ms', ...
    scale_ranges(1,1), scale_ranges(1,2), scale_ranges(2,1), scale_ranges(2,2));
end

function idx = match_numeric_axis(source, target)
idx = nan(size(source));
for i = 1:numel(source)
    [d, imin] = min(abs(target - source(i)));
    if d > 1e-6
        error('Could not match numeric axis value %.6f.', source(i));
    end
    idx(i) = imin;
end
end

function path = compute_path_length(traj, time, timerange)
time_idx = time >= timerange(1) & time <= timerange(2);
x = traj(:, time_idx, 1);
y = traj(:, time_idx, 2);
dx = diff(x, 1, 2);
dy = diff(y, 1, 2);
path = squeeze(sum(sqrt(dx.^2 + dy.^2), 2, 'omitnan'));
end

function amp = compute_lv3_reset_amplitude(traj, time, baseline_win, reset_win)
base_idx = time >= baseline_win(1) & time <= baseline_win(2);
reset_idx = time >= reset_win(1) & time <= reset_win(2);
lv3 = traj(:,:,2);
base = mean(lv3(:, base_idx), 2, 'omitnan');
lv3_bc = lv3 - base;
amp = max(abs(lv3_bc(:, reset_idx)), [], 2, 'omitnan');
end

function peak_speed = compute_peak_state_speed(traj, time, reset_win)
dt = median(diff(time));
tmid = time(1:end-1) + diff(time) ./ 2;
reset_idx = tmid >= reset_win(1) & tmid <= reset_win(2);
dx = diff(traj(:,:,1), 1, 2);
dy = diff(traj(:,:,2), 1, 2);
speed = sqrt(dx.^2 + dy.^2) ./ dt;
peak_speed = max(speed(:, reset_idx), [], 2, 'omitnan');
end

function outside_time = compute_time_outside_baseline_basin(traj, time, baseline_win, reset_win)
base_idx = time >= baseline_win(1) & time <= baseline_win(2);
reset_idx = time >= reset_win(1) & time <= reset_win(2);
dt = median(diff(time));
threshold_d2 = 5.991464547; % chi-square 95% threshold for 2-D state space
nSubj = size(traj, 1);
outside_time = nan(nSubj, 1);

for isub = 1:nSubj
    base_xy = squeeze(traj(isub, base_idx, 1:2));
    reset_xy = squeeze(traj(isub, reset_idx, 1:2));
    if size(base_xy, 1) < 3 || isempty(reset_xy)
        continue
    end

    center = mean(base_xy, 1, 'omitnan');
    C = cov(base_xy, 'omitrows');
    if any(isnan(C), 'all') || rcond(C) < 1e-8
        C = C + eye(2) * 1e-6;
    end

    dev = reset_xy - center;
    d2 = sum((dev / C) .* dev, 2);
    outside_time(isub) = sum(d2 > threshold_d2) * dt;
end
end

function [states, metrics] = compute_state_switch(traj, time, windows)
nSubj = size(traj, 1);
nWindows = size(windows, 1);
states = nan(nSubj, nWindows, 2);
for iwin = 1:nWindows
    idx = time >= windows(iwin, 1) & time <= windows(iwin, 2);
    if ~any(idx)
        error('Empty state-switch window %.2f to %.2f.', windows(iwin, 1), windows(iwin, 2));
    end
    states(:, iwin, :) = mean(traj(:, idx, :), 2, 'omitnan');
end

pre_state = reshape(states(:, 1, :), nSubj, 2);
quench_state = reshape(states(:, 2, :), nSubj, 2);
recovery_state = reshape(states(:, 3, :), nSubj, 2);
switch_vec = quench_state - pre_state;

metrics.lv1_switch = switch_vec(:, 1);
metrics.lv3_switch = switch_vec(:, 2);
metrics.switch_magnitude = sqrt(sum(switch_vec .^ 2, 2));
metrics.recovery_distance = sqrt(sum((recovery_state - pre_state) .^ 2, 2));
end

function [tbl, centroids] = compute_centroid_permutation(trajYA, trajOA, time, timerange, n_perm)
time_idx = time >= timerange(1) & time <= timerange(2);
if ~any(time_idx)
    error('Empty centroid time window %.2f to %.2f.', timerange(1), timerange(2));
end

centYA = [ ...
    mean(trajYA(:, time_idx, 1), 2, 'omitnan'), ...
    mean(trajYA(:, time_idx, 2), 2, 'omitnan')];
centOA = [ ...
    mean(trajOA(:, time_idx, 1), 2, 'omitnan'), ...
    mean(trajOA(:, time_idx, 2), 2, 'omitnan')];

centYA = centYA(all(~isnan(centYA), 2), :);
centOA = centOA(all(~isnan(centOA), 2), :);
nYA = size(centYA, 1);
nOA = size(centOA, 1);
if nYA < 2 || nOA < 2
    error('Need at least two valid subjects per group for centroid permutation.');
end

meanYA = mean(centYA, 1, 'omitnan');
meanOA = mean(centOA, 1, 'omitnan');
delta = meanOA - meanYA;
obs_dist = sqrt(sum(delta .^ 2));
if obs_dist > 0
    sep_axis = delta ./ obs_dist;
else
    sep_axis = [1 0];
end

all_cent = [centYA; centOA];
nTotal = size(all_cent, 1);
perm_dist = nan(n_perm, 1);
perm_lv1 = nan(n_perm, 1);
perm_lv3 = nan(n_perm, 1);

rng(11)
for iperm = 1:n_perm
    perm_idx = randperm(nTotal);
    idxYA = perm_idx(1:nYA);
    idxOA = perm_idx((nYA + 1):end);
    perm_delta = mean(all_cent(idxOA, :), 1, 'omitnan') - ...
        mean(all_cent(idxYA, :), 1, 'omitnan');
    perm_dist(iperm) = sqrt(sum(perm_delta .^ 2));
    perm_lv1(iperm) = abs(perm_delta(1));
    perm_lv3(iperm) = abs(perm_delta(2));
end

p_dist = (sum(perm_dist >= obs_dist) + 1) ./ (n_perm + 1);
p_lv1 = (sum(perm_lv1 >= abs(delta(1))) + 1) ./ (n_perm + 1);
p_lv3 = (sum(perm_lv3 >= abs(delta(2))) + 1) ./ (n_perm + 1);

tbl = table(nYA, nOA, timerange(1), timerange(2), ...
    meanYA(1), meanYA(2), meanOA(1), meanOA(2), ...
    delta(1), delta(2), obs_dist, p_dist, p_lv1, p_lv3, n_perm, ...
    'VariableNames', {'N_YA','N_OA','TimeStart','TimeEnd', ...
    'YA_LV1','YA_LV3','OA_LV1','OA_LV3', ...
    'Delta_LV1_OAminusYA','Delta_LV3_OAminusYA', ...
    'CentroidDistance','P_CentroidDistance','P_LV1','P_LV3','N_Permutations'});

centroids.YA = centYA;
centroids.OA = centOA;
centroids.separation_axis = sep_axis;
centroids.score_YA = centYA * sep_axis(:);
centroids.score_OA = centOA * sep_axis(:);
end

function [tbl, area_data] = compute_occupation_area_permutation(trajYA, trajOA, time, timerange, n_perm)
time_idx = time >= timerange(1) & time <= timerange(2);
if ~any(time_idx)
    error('Empty occupation time window %.2f to %.2f.', timerange(1), timerange(2));
end

trajYA = trajYA(:, time_idx, :);
trajOA = trajOA(:, time_idx, :);
nYA = size(trajYA, 1);
nOA = size(trajOA, 1);
if nYA < 2 || nOA < 2
    error('Need at least two valid subjects per group for occupation-area permutation.');
end

pooled_YA = occupation_area(trajYA);
pooled_OA = occupation_area(trajOA);
pooled_delta = pooled_OA - pooled_YA;
pooled_ratio = pooled_OA ./ pooled_YA;

subj_YA = subject_occupation_area(trajYA);
subj_OA = subject_occupation_area(trajOA);
subj_delta = mean(subj_OA, 'omitnan') - mean(subj_YA, 'omitnan');
subj_ratio = mean(subj_OA, 'omitnan') ./ mean(subj_YA, 'omitnan');

all_traj = cat(1, trajYA, trajOA);
all_subj = [subj_YA; subj_OA];
nTotal = size(all_traj, 1);
perm_pooled_delta = nan(n_perm, 1);
perm_subj_delta = nan(n_perm, 1);

rng(17)
for iperm = 1:n_perm
    perm_idx = randperm(nTotal);
    idxYA = perm_idx(1:nYA);
    idxOA = perm_idx((nYA + 1):end);

    perm_pooled_delta(iperm) = occupation_area(all_traj(idxOA,:,:)) - ...
        occupation_area(all_traj(idxYA,:,:));
    perm_subj_delta(iperm) = mean(all_subj(idxOA), 'omitnan') - ...
        mean(all_subj(idxYA), 'omitnan');
end

p_pooled_oa_gt_ya = (sum(perm_pooled_delta >= pooled_delta) + 1) ./ (n_perm + 1);
p_pooled_two_sided = (sum(abs(perm_pooled_delta) >= abs(pooled_delta)) + 1) ./ (n_perm + 1);
p_subj_oa_gt_ya = (sum(perm_subj_delta >= subj_delta) + 1) ./ (n_perm + 1);
p_subj_two_sided = (sum(abs(perm_subj_delta) >= abs(subj_delta)) + 1) ./ (n_perm + 1);

tbl = table(nYA, nOA, timerange(1), timerange(2), ...
    pooled_YA, pooled_OA, pooled_delta, pooled_ratio, ...
    p_pooled_oa_gt_ya, p_pooled_two_sided, ...
    mean(subj_YA, 'omitnan'), std(subj_YA, 0, 'omitnan') ./ sqrt(sum(~isnan(subj_YA))), ...
    mean(subj_OA, 'omitnan'), std(subj_OA, 0, 'omitnan') ./ sqrt(sum(~isnan(subj_OA))), ...
    subj_delta, subj_ratio, p_subj_oa_gt_ya, p_subj_two_sided, n_perm, ...
    'VariableNames', {'N_YA','N_OA','TimeStart','TimeEnd', ...
    'PooledArea_YA','PooledArea_OA','PooledArea_Delta_OAminusYA','PooledArea_Ratio_OAoverYA', ...
    'P_PooledArea_OAgtYA','P_PooledArea_TwoSided', ...
    'SubjectArea_YA_Mean','SubjectArea_YA_SEM','SubjectArea_OA_Mean','SubjectArea_OA_SEM', ...
    'SubjectArea_Delta_OAminusYA','SubjectArea_Ratio_OAoverYA', ...
    'P_SubjectArea_OAgtYA','P_SubjectArea_TwoSided','N_Permutations'});

area_data.subject_YA = subj_YA;
area_data.subject_OA = subj_OA;
end

function areas = subject_occupation_area(traj)
nSubj = size(traj, 1);
areas = nan(nSubj, 1);
for isub = 1:nSubj
    areas(isub) = occupation_area(traj(isub,:,:));
end
end

function area = occupation_area(traj)
xy = squeeze_trajectory_points(traj);
if size(xy, 1) < 3
    area = nan;
    return
end
C = cov(xy, 'omitrows');
if any(isnan(C), 'all')
    area = nan;
    return
end
area = pi * sqrt(max(det(C), 0));
end

function xy = squeeze_trajectory_points(traj)
if ndims(traj) ~= 3 || size(traj, 3) ~= 2
    error('Expected trajectory input as subject x time x 2.');
end
x = traj(:,:,1);
y = traj(:,:,2);
xy = [x(:), y(:)];
xy = xy(all(~isnan(xy), 2), :);
end

function tbl = compute_state_space_entropy_permutation(trajYA, trajOA, time, timerange, bin_counts, n_perm)
time_idx = time >= timerange(1) & time <= timerange(2);
if ~any(time_idx)
    error('Empty entropy time window %.2f to %.2f.', timerange(1), timerange(2));
end

trajYA = trajYA(:, time_idx, :);
trajOA = trajOA(:, time_idx, :);
nYA = size(trajYA, 1);
nOA = size(trajOA, 1);
if nYA < 2 || nOA < 2
    error('Need at least two valid subjects per group for state-space entropy permutation.');
end

all_xy = squeeze_trajectory_points(cat(1, trajYA, trajOA));
x_lim = padded_edges_limits(all_xy(:,1));
y_lim = padded_edges_limits(all_xy(:,2));

nRows = numel(bin_counts);
BinsPerAxis = bin_counts(:);
N_YA = repmat(nYA, nRows, 1);
N_OA = repmat(nOA, nRows, 1);
TimeStart = repmat(timerange(1), nRows, 1);
TimeEnd = repmat(timerange(2), nRows, 1);
N_Timepoints = repmat(nnz(time_idx), nRows, 1);
YA_Mean = nan(nRows, 1);
YA_SEM = nan(nRows, 1);
OA_Mean = nan(nRows, 1);
OA_SEM = nan(nRows, 1);
Delta_OAminusYA = nan(nRows, 1);
P_OAgtYA = nan(nRows, 1);
P_TwoSided = nan(nRows, 1);
P_Ttest = nan(nRows, 1);
T_Ttest = nan(nRows, 1);
EffectiveStates_YA = nan(nRows, 1);
EffectiveStates_OA = nan(nRows, 1);

rng(23)
for irow = 1:nRows
    n_bins = BinsPerAxis(irow);
    x_edges = linspace(x_lim(1), x_lim(2), n_bins + 1);
    y_edges = linspace(y_lim(1), y_lim(2), n_bins + 1);

    [entropy_YA, eff_YA] = subject_state_space_entropy(trajYA, x_edges, y_edges);
    [entropy_OA, eff_OA] = subject_state_space_entropy(trajOA, x_edges, y_edges);

    YA_Mean(irow) = mean(entropy_YA, 'omitnan');
    YA_SEM(irow) = std(entropy_YA, 0, 'omitnan') ./ sqrt(sum(~isnan(entropy_YA)));
    OA_Mean(irow) = mean(entropy_OA, 'omitnan');
    OA_SEM(irow) = std(entropy_OA, 0, 'omitnan') ./ sqrt(sum(~isnan(entropy_OA)));
    Delta_OAminusYA(irow) = OA_Mean(irow) - YA_Mean(irow);
    EffectiveStates_YA(irow) = mean(eff_YA, 'omitnan');
    EffectiveStates_OA(irow) = mean(eff_OA, 'omitnan');

    all_entropy = [entropy_YA; entropy_OA];
    perm_delta = nan(n_perm, 1);
    for iperm = 1:n_perm
        perm_idx = randperm(nYA + nOA);
        idxYA = perm_idx(1:nYA);
        idxOA = perm_idx((nYA + 1):end);
        perm_delta(iperm) = mean(all_entropy(idxOA), 'omitnan') - ...
            mean(all_entropy(idxYA), 'omitnan');
    end

    P_OAgtYA(irow) = (sum(perm_delta >= Delta_OAminusYA(irow)) + 1) ./ (n_perm + 1);
    P_TwoSided(irow) = (sum(abs(perm_delta) >= abs(Delta_OAminusYA(irow))) + 1) ./ (n_perm + 1);
    [~, p_t, ~, st] = ttest2(entropy_YA, entropy_OA, 'Vartype', 'unequal');
    P_Ttest(irow) = p_t;
    T_Ttest(irow) = st.tstat;
end

tbl = table(BinsPerAxis, N_YA, N_OA, TimeStart, TimeEnd, N_Timepoints, ...
    YA_Mean, YA_SEM, OA_Mean, OA_SEM, Delta_OAminusYA, ...
    P_OAgtYA, P_TwoSided, P_Ttest, T_Ttest, ...
    EffectiveStates_YA, EffectiveStates_OA, repmat(n_perm, nRows, 1), ...
    'VariableNames', {'BinsPerAxis','N_YA','N_OA','TimeStart','TimeEnd','N_Timepoints', ...
    'Entropy_YA_Mean','Entropy_YA_SEM','Entropy_OA_Mean','Entropy_OA_SEM', ...
    'Entropy_Delta_OAminusYA','P_Entropy_OAgtYA','P_Entropy_TwoSided', ...
    'P_Entropy_Ttest','T_Entropy_Ttest', ...
    'EffectiveStates_YA_Mean','EffectiveStates_OA_Mean','N_Permutations'});
end

function lim = padded_edges_limits(vals)
vals = vals(~isnan(vals));
lo = min(vals);
hi = max(vals);
if lo == hi
    pad = max(0.01, abs(lo) * 0.05);
else
    pad = max(0.01, (hi - lo) * 0.001);
end
lim = [lo - pad, hi + pad];
end

function [entropy_vals, effective_states] = subject_state_space_entropy(traj, x_edges, y_edges)
nSubj = size(traj, 1);
n_bins = numel(x_edges) - 1;
entropy_vals = nan(nSubj, 1);
effective_states = nan(nSubj, 1);

for isub = 1:nSubj
    xy = squeeze_trajectory_points(traj(isub,:,:));
    if isempty(xy)
        continue
    end
    xbin = discretize(xy(:,1), x_edges);
    ybin = discretize(xy(:,2), y_edges);
    valid = ~isnan(xbin) & ~isnan(ybin);
    if ~any(valid)
        continue
    end
    bin_id = xbin(valid) + (ybin(valid) - 1) .* n_bins;
    counts = accumarray(bin_id, 1, [n_bins * n_bins 1]);
    p = counts ./ sum(counts);
    p = p(p > 0);
    h = -sum(p .* log2(p));
    entropy_vals(isub) = h ./ log2(n_bins * n_bins);
    effective_states(isub) = 2 .^ h;
end
end

function tbl = compute_state_transition_entropy_permutation(trajYA, trajOA, time, timerange, bin_counts, n_perm)
time_idx = time >= timerange(1) & time <= timerange(2);
if ~any(time_idx)
    error('Empty transition-entropy time window %.2f to %.2f.', timerange(1), timerange(2));
end

trajYA = trajYA(:, time_idx, :);
trajOA = trajOA(:, time_idx, :);
nYA = size(trajYA, 1);
nOA = size(trajOA, 1);
if nYA < 2 || nOA < 2
    error('Need at least two valid subjects per group for transition entropy permutation.');
end

all_xy = squeeze_trajectory_points(cat(1, trajYA, trajOA));
x_lim = padded_edges_limits(all_xy(:,1));
y_lim = padded_edges_limits(all_xy(:,2));

nRows = numel(bin_counts);
BinsPerAxis = bin_counts(:);
N_YA = repmat(nYA, nRows, 1);
N_OA = repmat(nOA, nRows, 1);
TimeStart = repmat(timerange(1), nRows, 1);
TimeEnd = repmat(timerange(2), nRows, 1);
N_Transitions = repmat(max(nnz(time_idx) - 1, 0), nRows, 1);
YA_Mean = nan(nRows, 1);
YA_SEM = nan(nRows, 1);
OA_Mean = nan(nRows, 1);
OA_SEM = nan(nRows, 1);
Delta_OAminusYA = nan(nRows, 1);
P_OAgtYA = nan(nRows, 1);
P_TwoSided = nan(nRows, 1);
P_Ttest = nan(nRows, 1);
T_Ttest = nan(nRows, 1);
EffectiveTransitions_YA = nan(nRows, 1);
EffectiveTransitions_OA = nan(nRows, 1);

rng(29)
for irow = 1:nRows
    n_bins = BinsPerAxis(irow);
    x_edges = linspace(x_lim(1), x_lim(2), n_bins + 1);
    y_edges = linspace(y_lim(1), y_lim(2), n_bins + 1);

    [entropy_YA, eff_YA] = subject_transition_entropy(trajYA, x_edges, y_edges);
    [entropy_OA, eff_OA] = subject_transition_entropy(trajOA, x_edges, y_edges);

    YA_Mean(irow) = mean(entropy_YA, 'omitnan');
    YA_SEM(irow) = std(entropy_YA, 0, 'omitnan') ./ sqrt(sum(~isnan(entropy_YA)));
    OA_Mean(irow) = mean(entropy_OA, 'omitnan');
    OA_SEM(irow) = std(entropy_OA, 0, 'omitnan') ./ sqrt(sum(~isnan(entropy_OA)));
    Delta_OAminusYA(irow) = OA_Mean(irow) - YA_Mean(irow);
    EffectiveTransitions_YA(irow) = mean(eff_YA, 'omitnan');
    EffectiveTransitions_OA(irow) = mean(eff_OA, 'omitnan');

    all_entropy = [entropy_YA; entropy_OA];
    perm_delta = nan(n_perm, 1);
    for iperm = 1:n_perm
        perm_idx = randperm(nYA + nOA);
        idxYA = perm_idx(1:nYA);
        idxOA = perm_idx((nYA + 1):end);
        perm_delta(iperm) = mean(all_entropy(idxOA), 'omitnan') - ...
            mean(all_entropy(idxYA), 'omitnan');
    end

    P_OAgtYA(irow) = (sum(perm_delta >= Delta_OAminusYA(irow)) + 1) ./ (n_perm + 1);
    P_TwoSided(irow) = (sum(abs(perm_delta) >= abs(Delta_OAminusYA(irow))) + 1) ./ (n_perm + 1);
    [~, p_t, ~, st] = ttest2(entropy_YA, entropy_OA, 'Vartype', 'unequal');
    P_Ttest(irow) = p_t;
    T_Ttest(irow) = st.tstat;
end

tbl = table(BinsPerAxis, N_YA, N_OA, TimeStart, TimeEnd, N_Transitions, ...
    YA_Mean, YA_SEM, OA_Mean, OA_SEM, Delta_OAminusYA, ...
    P_OAgtYA, P_TwoSided, P_Ttest, T_Ttest, ...
    EffectiveTransitions_YA, EffectiveTransitions_OA, repmat(n_perm, nRows, 1), ...
    'VariableNames', {'BinsPerAxis','N_YA','N_OA','TimeStart','TimeEnd','N_Transitions', ...
    'TransitionEntropy_YA_Mean','TransitionEntropy_YA_SEM', ...
    'TransitionEntropy_OA_Mean','TransitionEntropy_OA_SEM', ...
    'TransitionEntropy_Delta_OAminusYA','P_TransitionEntropy_OAgtYA', ...
    'P_TransitionEntropy_TwoSided','P_TransitionEntropy_Ttest','T_TransitionEntropy_Ttest', ...
    'EffectiveTransitions_YA_Mean','EffectiveTransitions_OA_Mean','N_Permutations'});
end

function [entropy_vals, effective_transitions] = subject_transition_entropy(traj, x_edges, y_edges)
nSubj = size(traj, 1);
n_bins = numel(x_edges) - 1;
n_states = n_bins * n_bins;
entropy_vals = nan(nSubj, 1);
effective_transitions = nan(nSubj, 1);

for isub = 1:nSubj
    xy = squeeze_trajectory_points(traj(isub,:,:));
    if size(xy, 1) < 2
        continue
    end
    xbin = discretize(xy(:,1), x_edges);
    ybin = discretize(xy(:,2), y_edges);
    valid = ~isnan(xbin) & ~isnan(ybin);
    state = nan(size(xbin));
    state(valid) = xbin(valid) + (ybin(valid) - 1) .* n_bins;

    from_state = state(1:end-1);
    to_state = state(2:end);
    valid_trans = ~isnan(from_state) & ~isnan(to_state);
    if ~any(valid_trans)
        continue
    end

    trans_id = from_state(valid_trans) + (to_state(valid_trans) - 1) .* n_states;
    counts = accumarray(trans_id, 1, [n_states * n_states 1]);
    p = counts ./ sum(counts);
    p = p(p > 0);
    h = -sum(p .* log2(p));
    entropy_vals(isub) = h ./ log2(n_states * n_states);
    effective_transitions(isub) = 2 .^ h;
end
end

function [tbl, metric_data] = compute_escape_return_permutation( ...
    trajYA, trajOA, time, baseline_win, escape_win, post_win, n_perm)
metrics_YA = subject_escape_return_metrics(trajYA, time, baseline_win, escape_win, post_win);
metrics_OA = subject_escape_return_metrics(trajOA, time, baseline_win, escape_win, post_win);

metric_names = ["Peak basin distance"; "Time outside basin"; "Return latency"; "Max post-blink speed"];
directions = ["OA>YA"; "OA>YA"; "OA>YA"; "OA>YA"];

tbl = compare_group_metric_matrix(metrics_YA, metrics_OA, metric_names, directions, n_perm);
tbl.BaselineStart = repmat(baseline_win(1), height(tbl), 1);
tbl.BaselineEnd = repmat(baseline_win(2), height(tbl), 1);
tbl.EscapeStart = repmat(escape_win(1), height(tbl), 1);
tbl.EscapeEnd = repmat(escape_win(2), height(tbl), 1);
tbl.PostStart = repmat(post_win(1), height(tbl), 1);
tbl.PostEnd = repmat(post_win(2), height(tbl), 1);

metric_data.YA = metrics_YA;
metric_data.OA = metrics_OA;
metric_data.metric_names = metric_names;
end

function metrics = subject_escape_return_metrics(traj, time, baseline_win, escape_win, post_win)
base_idx = time >= baseline_win(1) & time <= baseline_win(2);
escape_idx = time >= escape_win(1) & time <= escape_win(2);
post_idx = time >= post_win(1) & time <= post_win(2);
dt = median(diff(time));
threshold_d2 = 5.991464547; % chi-square 95% threshold for 2-D state space
nSubj = size(traj, 1);
metrics = nan(nSubj, 4);

for isub = 1:nSubj
    xy = squeeze(traj(isub,:,1:2));
    if size(xy, 2) ~= 2
        xy = squeeze(xy)';
    end
    base_xy = xy(base_idx, :);
    if size(base_xy, 1) < 3
        continue
    end

    center = mean(base_xy, 1, 'omitnan');
    C = cov(base_xy, 'omitrows');
    if any(isnan(C), 'all') || rcond(C) < 1e-8
        C = C + eye(2) * 1e-6;
    end

    dev = xy - center;
    d2 = sum((dev / C) .* dev, 2);
    dist = sqrt(max(d2, 0));

    metrics(isub, 1) = max(dist(escape_idx), [], 'omitnan');
    metrics(isub, 2) = sum(d2(post_idx) > threshold_d2) .* dt;
    metrics(isub, 3) = compute_return_latency(time, d2, post_idx, threshold_d2, dt);
    metrics(isub, 4) = compute_subject_peak_speed(traj(isub,:,:), time, escape_win);
end
end

function latency = compute_return_latency(time, d2, post_idx, threshold_d2, dt)
post_times = time(post_idx);
post_d2 = d2(post_idx);
outside = post_d2 > threshold_d2;
if isempty(post_d2) || ~any(outside)
    latency = 0;
    return
end

first_out = find(outside, 1, 'first');
inside_after = find(post_d2(first_out:end) <= threshold_d2, 1, 'first');
if isempty(inside_after)
    latency = (post_times(end) - post_times(first_out)) + dt;
else
    return_idx = first_out + inside_after - 1;
    latency = post_times(return_idx) - post_times(first_out);
end
end

function peak_speed = compute_subject_peak_speed(traj, time, win)
dt = median(diff(time));
tmid = time(1:end-1) + diff(time) ./ 2;
idx = tmid >= win(1) & tmid <= win(2);
xy = squeeze(traj);
if size(xy, 2) ~= 2
    xy = xy';
end
dxy = diff(xy, 1, 1);
speed = sqrt(sum(dxy .^ 2, 2)) ./ dt;
peak_speed = max(speed(idx), [], 'omitnan');
end

function tbl = compare_group_metric_matrix(valsYA, valsOA, metric_names, directions, n_perm)
nMetrics = numel(metric_names);
Metric = metric_names(:);
Direction = directions(:);
N_YA = nan(nMetrics, 1);
N_OA = nan(nMetrics, 1);
YA_Mean = nan(nMetrics, 1);
YA_SEM = nan(nMetrics, 1);
OA_Mean = nan(nMetrics, 1);
OA_SEM = nan(nMetrics, 1);
Delta_OAminusYA = nan(nMetrics, 1);
P_OAgtYA = nan(nMetrics, 1);
P_TwoSided = nan(nMetrics, 1);
P_Ttest = nan(nMetrics, 1);
T_Ttest = nan(nMetrics, 1);
CohenD_OAminusYA = nan(nMetrics, 1);

rng(31)
for imetric = 1:nMetrics
    y = valsYA(:, imetric);
    o = valsOA(:, imetric);
    valid_y = ~isnan(y);
    valid_o = ~isnan(o);
    y = y(valid_y);
    o = o(valid_o);
    N_YA(imetric) = numel(y);
    N_OA(imetric) = numel(o);
    YA_Mean(imetric) = mean(y, 'omitnan');
    YA_SEM(imetric) = std(y, 0, 'omitnan') ./ sqrt(numel(y));
    OA_Mean(imetric) = mean(o, 'omitnan');
    OA_SEM(imetric) = std(o, 0, 'omitnan') ./ sqrt(numel(o));
    Delta_OAminusYA(imetric) = OA_Mean(imetric) - YA_Mean(imetric);

    all_vals = [y; o];
    nY = numel(y);
    nO = numel(o);
    perm_delta = nan(n_perm, 1);
    for iperm = 1:n_perm
        perm_idx = randperm(nY + nO);
        idxY = perm_idx(1:nY);
        idxO = perm_idx((nY + 1):end);
        perm_delta(iperm) = mean(all_vals(idxO), 'omitnan') - ...
            mean(all_vals(idxY), 'omitnan');
    end

    P_OAgtYA(imetric) = (sum(perm_delta >= Delta_OAminusYA(imetric)) + 1) ./ (n_perm + 1);
    P_TwoSided(imetric) = (sum(abs(perm_delta) >= abs(Delta_OAminusYA(imetric))) + 1) ./ (n_perm + 1);
    [~, p_t, ~, st] = ttest2(y, o, 'Vartype', 'unequal');
    P_Ttest(imetric) = p_t;
    T_Ttest(imetric) = st.tstat;

    sd_pooled = sqrt(((nY - 1) * var(y, 'omitnan') + ...
        (nO - 1) * var(o, 'omitnan')) ./ (nY + nO - 2));
    CohenD_OAminusYA(imetric) = Delta_OAminusYA(imetric) ./ sd_pooled;
end

tbl = table(Metric, Direction, N_YA, N_OA, YA_Mean, YA_SEM, OA_Mean, OA_SEM, ...
    Delta_OAminusYA, P_OAgtYA, P_TwoSided, P_Ttest, T_Ttest, CohenD_OAminusYA, ...
    repmat(n_perm, nMetrics, 1), ...
    'VariableNames', {'Metric','Direction','N_YA','N_OA','YA_Mean','YA_SEM', ...
    'OA_Mean','OA_SEM','Delta_OAminusYA','P_OAgtYA','P_TwoSided', ...
    'P_Ttest','T_Ttest','CohenD_OAminusYA','N_Permutations'});
end

function out = window_mean(tc, time, win)
idx = time >= win(1) & time <= win(2);
if ~any(idx)
    error('Empty time window %.2f to %.2f.', win(1), win(2));
end
out = mean(tc(:, idx), 2, 'omitnan');
end

function out = window_min(tc, time, win)
idx = time >= win(1) & time <= win(2);
if ~any(idx)
    error('Empty time window %.2f to %.2f.', win(1), win(2));
end
out = min(tc(:, idx), [], 2, 'omitnan');
end

function yl = mean_sem_ylim(groups, pad_frac)
vals = [];
for igroup = 1:numel(groups)
    tc = groups{igroup};
    m = mean(tc, 1, 'omitnan');
    sem = std(tc, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(tc), 1));
    vals = [vals, m - sem, m + sem]; %#ok<AGROW>
end
yl = padded_range(vals, pad_frac);
end

function yl = padded_range(vals, pad_frac)
vals = vals(~isnan(vals));
if isempty(vals)
    yl = [-1 1];
    return
end
lo = min(vals);
hi = max(vals);
if lo == hi
    pad = max(0.01, abs(lo) * 0.05);
else
    pad = max(0.01, (hi - lo) * pad_frac);
end
yl = [lo - pad, hi + pad];
end

function tc_bc = baseline_correct(tc, time, baseline_win)
base = window_mean(tc, time, baseline_win);
tc_bc = tc - base;
end

function plot_entropy_timecourse(ax, time, tcYA, tcOA, colYA, colOA, ttl, ylab, rise_base, rise_pre, quench_pre, quench_post)
hold(ax, 'on')
hYA = plot_shaded(ax, time, tcYA, colYA);
hOA = plot_shaded(ax, time, tcOA, colOA);

mYA = mean(tcYA, 1, 'omitnan');
mOA = mean(tcOA, 1, 'omitnan');
semYA = std(tcYA, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(tcYA), 1));
semOA = std(tcOA, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(tcOA), 1));
yl = [min([mYA - semYA, mOA - semOA]) max([mYA + semYA, mOA + semOA])];
ypad = max(0.01, diff(yl) * 0.45);
ylim(ax, [yl(1) - ypad yl(2) + ypad]);

shade_time_window(ax, rise_base, [0.92 0.92 0.92]);
shade_time_window(ax, rise_pre, [0.96 0.90 0.78]);
shade_time_window(ax, quench_pre, [0.96 0.90 0.78]);
shade_time_window(ax, quench_post, [0.84 0.90 0.98]);
xline(ax, 0, 'k:', 'LineWidth', 0.9);

uistack(hYA, 'top');
uistack(hOA, 'top');
xlim(ax, [-1.00 1.00]);
xlabel(ax, 'Time from blink peak (s)');
ylabel(ax, ylab);
title(ax, ttl, 'FontWeight', 'bold');
text(ax, 0.93, mYA(end), 'YA', 'Color', colYA, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
text(ax, 0.93, mOA(end), 'OA', 'Color', colOA, 'FontWeight', 'bold', ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');
format_axes(ax);
end

function shade_time_window(ax, win, color)
yl = ax.YLim;
p = patch(ax, [win(1) win(2) win(2) win(1)], [yl(1) yl(1) yl(2) yl(2)], color, ...
    'EdgeColor', 'none', 'FaceAlpha', 0.32, 'HandleVisibility', 'off');
uistack(p, 'bottom');
end

function h = plot_shaded(ax, time, tc, col)
m = mean(tc, 1, 'omitnan');
sem = std(tc, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(tc), 1));
fill(ax, [time fliplr(time)], [m + sem fliplr(m - sem)], col, ...
    'FaceAlpha', 0.18, 'EdgeColor', 'none', 'HandleVisibility', 'off');
h = plot(ax, time, m, 'Color', col, 'LineWidth', 1.8);
end

function plot_state_switch(ax, statesYA, statesOA, colYA, colOA)
hold(ax, 'on')
[meansYA, hYA] = plot_switch_group(ax, statesYA, colYA);
[meansOA, hOA] = plot_switch_group(ax, statesOA, colOA);

xline(ax, 0, '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.5, 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel(ax, 'LV1 expression');
ylabel(ax, 'LV3 expression');

hPre = plot(ax, nan, nan, 'o', 'Color', 'k', 'MarkerFaceColor', 'w', ...
    'MarkerSize', 4.0, 'LineWidth', 0.7);
hQuench = plot(ax, nan, nan, 'd', 'Color', 'k', 'MarkerFaceColor', 'w', ...
    'MarkerSize', 4.0, 'LineWidth', 0.7);
hRecovery = plot(ax, nan, nan, 's', 'Color', 'k', 'MarkerFaceColor', 'w', ...
    'MarkerSize', 4.0, 'LineWidth', 0.7);

xvals = [meansYA(:,1); meansOA(:,1)];
yvals = [meansYA(:,2); meansOA(:,2)];
xlim(ax, padded_range(xvals(:), 0.10));
ylim(ax, padded_range(yvals(:), 0.10));

format_axes(ax);
lgd = legend(ax, [hYA hOA hPre hQuench hRecovery], ...
    {'YA','OA','pre','quench','recovery'}, 'Location', 'northwest', ...
    'Box', 'off', 'FontSize', 5.5);
lgd.ItemTokenSize = [8 5];
end

function [means, hLine] = plot_switch_group(ax, states, col)
nWindows = size(states, 2);
means = nan(nWindows, 2);
sems = nan(nWindows, 2);
for iwin = 1:nWindows
    pts = reshape(states(:, iwin, :), [], 2);
    means(iwin,:) = mean(pts, 1, 'omitnan');
    sems(iwin,:) = std(pts, 0, 1, 'omitnan') ./ sqrt(sum(all(~isnan(pts), 2)));
end

light_col = lighten_color(col, 0.62);
for iwin = 1:nWindows
    draw_sem_cross(ax, means(iwin,:), sems(iwin,:), light_col);
end

hLine = plot(ax, means(:,1), means(:,2), '-', 'Color', col, 'LineWidth', 2.4, ...
    'HandleVisibility', 'off');
for iseg = 1:(nWindows - 1)
    delta = means(iseg + 1,:) - means(iseg,:);
    quiver(ax, means(iseg,1), means(iseg,2), delta(1), delta(2), 0, ...
        'Color', col, 'LineWidth', 0.85, 'MaxHeadSize', 0.65, ...
        'AutoScale', 'off', 'HandleVisibility', 'off');
end

markers = {'o', 'd', 's'};
for iwin = 1:nWindows
    plot(ax, means(iwin,1), means(iwin,2), markers{iwin}, ...
        'Color', col, 'MarkerFaceColor', col, 'MarkerEdgeColor', 'k', ...
        'MarkerSize', 6.0, 'LineWidth', 0.7, 'HandleVisibility', 'off');
end
end

function draw_sem_cross(ax, mean_xy, sem_xy, col)
plot(ax, [mean_xy(1) - sem_xy(1), mean_xy(1) + sem_xy(1)], ...
    [mean_xy(2), mean_xy(2)], '-', 'Color', col, 'LineWidth', 0.75, ...
    'HandleVisibility', 'off');
plot(ax, [mean_xy(1), mean_xy(1)], ...
    [mean_xy(2) - sem_xy(2), mean_xy(2) + sem_xy(2)], '-', ...
    'Color', col, 'LineWidth', 0.75, 'HandleVisibility', 'off');
end

function plot_state_space(ax, trajYA, trajOA, time, colYA, colOA, timerange)
hold(ax, 'on')
time_idx = time >= timerange(1) & time <= timerange(2);
xYA = trajYA(:, time_idx, 1);
yYA = trajYA(:, time_idx, 2);
xOA = trajOA(:, time_idx, 1);
yOA = trajOA(:, time_idx, 2);
tplot = time(time_idx);

ellYA = plot_cov_ellipse(ax, xYA(:), yYA(:), colYA);
ellOA = plot_cov_ellipse(ax, xOA(:), yOA(:), colOA);

mxYA = mean(xYA, 1, 'omitnan');
myYA = mean(yYA, 1, 'omitnan');
mxOA = mean(xOA, 1, 'omitnan');
myOA = mean(yOA, 1, 'omitnan');

plot_mean_trajectory(ax, mxYA, myYA, tplot, colYA);
plot_mean_trajectory(ax, mxOA, myOA, tplot, colOA);

xline(ax, 0, '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.5, 'HandleVisibility', 'off');
yline(ax, 0, '-', 'Color', [0.82 0.82 0.82], 'LineWidth', 0.5, 'HandleVisibility', 'off');
xlabel(ax, 'LV1 expression');
ylabel(ax, 'LV3 expression');
ylim(ax, padded_range([ellYA(2,:), ellOA(2,:), myYA, myOA], 0.05));
format_axes(ax);
end

function plot_mean_trajectory(ax, x, y, time, col)
plot(ax, x, y, '-', 'Color', col, 'LineWidth', 2.6);
mark_idx = 1:5:numel(time);
plot(ax, x(mark_idx), y(mark_idx), 'o', 'Color', col, 'MarkerFaceColor', 'w', ...
    'MarkerSize', 3.8, 'LineWidth', 0.9, 'HandleVisibility', 'off');
[~, blink_idx] = min(abs(time));
plot(ax, x(1), y(1), '^', 'Color', col, 'MarkerFaceColor', col, ...
    'MarkerSize', 5.5, 'HandleVisibility', 'off');
plot(ax, x(blink_idx), y(blink_idx), 'o', 'Color', col, 'MarkerFaceColor', col, ...
    'MarkerEdgeColor', 'k', 'MarkerSize', 5.5, 'HandleVisibility', 'off');
plot(ax, x(end), y(end), 's', 'Color', col, 'MarkerFaceColor', col, ...
    'MarkerSize', 5.2, 'HandleVisibility', 'off');

arrow_idx = 4:8:(numel(time)-1);
quiver(ax, x(arrow_idx), y(arrow_idx), diff(x([arrow_idx; arrow_idx + 1])), ...
    diff(y([arrow_idx; arrow_idx + 1])), 0, 'Color', col, 'LineWidth', 0.8, ...
    'MaxHeadSize', 0.7, 'AutoScale', 'off', 'HandleVisibility', 'off');
end

function ellipse = plot_cov_ellipse(ax, x, y, col)
valid = ~isnan(x) & ~isnan(y);
x = x(valid);
y = y(valid);
if numel(x) < 3
    ellipse = nan(2, 0);
    return
end
C = cov([x y]);
[V, D] = eig(C);
theta = linspace(0, 2*pi, 160);
circle = [cos(theta); sin(theta)];
scale = sqrt(0.2107); % about 10% for a 2D normal cloud
ellipse = scale * V * sqrt(D) * circle;
ellipse(1,:) = ellipse(1,:) + mean(x);
ellipse(2,:) = ellipse(2,:) + mean(y);
fill(ax, ellipse(1,:), ellipse(2,:), col, 'FaceAlpha', 0.025, ...
    'EdgeColor', col, 'LineWidth', 0.55, 'HandleVisibility', 'off');
end

function plot_metric_panel(ax, YA, OA, colYA, colOA, ttl, ylab, pval)
hold(ax, 'on')
rng(1)
xYA = 1 + (rand(size(YA)) - 0.5) * 0.16;
xOA = 2 + (rand(size(OA)) - 0.5) * 0.16;
scatter(ax, xYA, YA, 15, colYA, 'filled', 'MarkerFaceAlpha', 0.52, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.25);
scatter(ax, xOA, OA, 15, colOA, 'filled', 'MarkerFaceAlpha', 0.52, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 0.25);

mYA = mean(YA, 'omitnan');
mOA = mean(OA, 'omitnan');
semYA = std(YA, 0, 'omitnan') ./ sqrt(sum(~isnan(YA)));
semOA = std(OA, 0, 'omitnan') ./ sqrt(sum(~isnan(OA)));
bar(ax, 1, mYA, 0.48, 'FaceColor', colYA, 'FaceAlpha', 0.28, 'EdgeColor', 'none');
bar(ax, 2, mOA, 0.48, 'FaceColor', colOA, 'FaceAlpha', 0.28, 'EdgeColor', 'none');
errorbar(ax, [1 2], [mYA mOA], [semYA semOA], 'k', 'LineStyle', 'none', ...
    'LineWidth', 0.9, 'CapSize', 5);
plot(ax, [1 2], [mYA mOA], 'k-', 'LineWidth', 0.7);

xlim(ax, [0.55 2.45]);
yl = padded_range([YA(:); OA(:); 0], 0.06);
if min([YA(:); OA(:)]) >= 0
    yl(1) = 0;
end
ylim(ax, yl);
if ~isnan(pval)
    txt_weight = 'normal';
    if pval < 0.05
        txt_weight = 'bold';
    end
    text(ax, 1.5, yl(2) - 0.10 * diff(yl), p_string(pval), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', 7, 'FontWeight', txt_weight);
end
xticks(ax, [1 2]);
xticklabels(ax, {'YA','OA'});
ylabel(ax, ylab);
title(ax, ttl, 'FontWeight', 'bold');
format_axes(ax);
end

function tbl = make_metric_table(metric_names, vals_YA, vals_OA)
n = numel(metric_names);
YA_Mean = nan(n,1);
YA_SEM = nan(n,1);
OA_Mean = nan(n,1);
OA_SEM = nan(n,1);
T_Group = nan(n,1);
P_Group = nan(n,1);
CohenD_OAminusYA = nan(n,1);
for i = 1:n
    y = vals_YA{i};
    o = vals_OA{i};
    YA_Mean(i) = mean(y, 'omitnan');
    YA_SEM(i) = std(y, 0, 'omitnan') ./ sqrt(sum(~isnan(y)));
    OA_Mean(i) = mean(o, 'omitnan');
    OA_SEM(i) = std(o, 0, 'omitnan') ./ sqrt(sum(~isnan(o)));
    [~, p, ~, st] = ttest2(y, o, 'Vartype', 'unequal');
    T_Group(i) = st.tstat;
    P_Group(i) = p;
    sd_pooled = sqrt(((sum(~isnan(y))-1) * var(y, 'omitnan') + ...
        (sum(~isnan(o))-1) * var(o, 'omitnan')) ./ ...
        (sum(~isnan(y)) + sum(~isnan(o)) - 2));
    CohenD_OAminusYA(i) = (OA_Mean(i) - YA_Mean(i)) ./ sd_pooled;
end
tbl = table(metric_names, YA_Mean, YA_SEM, OA_Mean, OA_SEM, ...
    T_Group, P_Group, CohenD_OAminusYA, ...
    'VariableNames', {'Metric','YA_Mean','YA_SEM','OA_Mean','OA_SEM', ...
    'T_Group','P_Group','CohenD_OAminusYA'});
end

function shade = lighten_color(col, amount)
shade = col + (1 - col) .* amount;
end

function txt = p_string(p)
if isnan(p)
    txt = 'p = n/a';
elseif p < .001
    txt = 'p < .001';
else
    txt = sprintf('p = %.3f', p);
end
end

function format_axes(ax)
box(ax, 'off')
set(ax, 'FontSize', 8, 'LineWidth', 0.8, 'TickDir', 'out', ...
    'XColor', 'k', 'YColor', 'k');
if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
    ax.Toolbar.Visible = 'off';
end
end
