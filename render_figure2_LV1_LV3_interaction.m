%% Figure 2: LV1, LV3, and their regime-dependent interaction
%
% Panels:
%   A) LV1 bootstrap-ratio topography and time-by-timescale map
%   B) LV1 brain-score/raw behavior correlations
%   C) LV3 bootstrap-ratio map with two transient topography windows
%   D) LV3 brain-score/raw behavior correlations
%   E) LV1 vs LV3 brain-score scatter
%   F) LV1 x LV3 behavior coefficient across cognitive regimes

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

folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

if ~exist(plotfolder, 'dir')
    mkdir(plotfolder);
end

if exist('stat_mse_2group_behav', 'var')
    fprintf('Using stat_mse_2group_behav already present in the MATLAB workspace.\n');
else
    stat_file = fullfile(plotfolder, 'behavPLS_2group_blink.mat');
    if ~exist(stat_file, 'file')
        stat_file = fullfile(folder, 'behavPLS_2group_blink.mat');
    end
    fprintf('Loading stat_mse_2group_behav from %s\n', char(stat_file));
    S = load(stat_file, 'stat_mse_2group_behav');
    stat_mse_2group_behav = S.stat_mse_2group_behav;
end
stat_mse_2group_behav = orient_behav_pls_lvs(stat_mse_2group_behav);

if exist('results', 'var') && istable(results) && ...
        all(ismember({'B_LV1xLV3_YA','SE_LV1xLV3_YA'}, results.Properties.VariableNames))
    interaction_results = results;
else
    results_file = fullfile(plotfolder, 'LV1xLV3_predict_behavior_results.mat');
    if ~exist(results_file, 'file')
        test_LV1xLV3_predict_behavior
    end
    Sres = load(results_file, 'results');
    interaction_results = Sres.results;
end

if ~exist('interaction_results', 'var')
    test_LV1xLV3_predict_behavior
    interaction_results = results;
end

cmap = cbrewer('div','RdBu',256);
cmap = flipud(cmap);
col_YA = [0.85 0.25 0.25];
col_OA = [0.20 0.45 0.85];
col_YA_light = [0.92 0.55 0.55];
col_OA_light = [0.52 0.68 0.92];

f = figure('Color', 'w');
f.Units = 'centimeters';
f.Position = [2 2 18 19];
colormap(f, cmap);

map_pos_1 = [0.075 0.705 0.335 0.255];
bar_pos_1 = [0.600 0.735 0.325 0.190];
map_pos_3 = [0.075 0.390 0.335 0.255];
bar_pos_3 = [0.600 0.420 0.325 0.190];
scatter_pos = [0.105 0.090 0.325 0.205];
interaction_pos = [0.610 0.090 0.325 0.205];

lv1_window.timerange = [-1.00 -0.25];
lv1_window.freqrange = [];
lv1_window.label = '-1.00 to -0.25 s';
lv3_windows = lv1_window;
lv3_windows(1).timerange = [-1.25 -0.30];
lv3_windows(1).freqrange = [];
lv3_windows(1).label = 'Pre';
lv3_windows(2).timerange = [-0.20 0.40];
lv3_windows(2).freqrange = [20 60];
lv3_windows(2).label = 'Peri';

plot_lv_map(f, map_pos_1, stat_mse_2group_behav, 1, lv1_window, cmap, 'A  LV1 feature reliability');
plot_behav_corr_panel(f, bar_pos_1, stat_mse_2group_behav, 1, col_YA_light, col_OA_light, col_YA, col_OA, 'B  LV1 brain score vs behavior', false);

plot_lv_map(f, map_pos_3, stat_mse_2group_behav, 3, lv3_windows, cmap, 'C  LV3 feature reliability');
plot_behav_corr_panel(f, bar_pos_3, stat_mse_2group_behav, 3, col_YA_light, col_OA_light, col_YA, col_OA, 'D  LV3 brain score vs behavior', false);

plot_lv1_lv3_scatter_panel(f, scatter_pos, stat_mse_2group_behav, col_YA, col_OA);
plot_interaction_panel(f, interaction_pos, interaction_results, col_YA, col_OA);

outfile = 'Figure2_LV1_LV3_interaction';
force_map_colormaps(f, cmap);
exportFigure(f, fullfile(plotfolder, [outfile '.pdf']), 'Resolution', 600, 'FontName', 'Arial');
force_map_colormaps(f, cmap);
exportFigure(f, fullfile(plotfolder, [outfile '.png']), 'Resolution', 600, 'FontName', 'Arial');

fprintf('Saved %s outputs to %s\n', outfile, char(plotfolder));

%% Local helpers

function force_map_colormaps(fig, cmap)
mapAxes = findall(fig, 'Type', 'axes');
for iax = 1:numel(mapAxes)
    hasImage = ~isempty(findall(mapAxes(iax), 'Type', 'image'));
    hasSurface = ~isempty(findall(mapAxes(iax), 'Type', 'surface'));
    if hasImage || hasSurface
        colormap(mapAxes(iax), cmap);
    end
end
colormap(fig, cmap);
end

function stat = orient_behav_pls_lvs(stat)
% PLS signs are arbitrary. Match the manuscript convention used in
% LV1_BehavPLS_domain.png: LV1 should run from negative young/stable-domain
% correlations toward older/flexible-domain compensation.
ord = [5 4 3 1 2];
for lv = 1
    behav_YA = stat.results.stacked_behavdata(1:20, ord);
    brainscore_YA = stat.brainscores{1}(:, lv);
    r_stable_YA = corr(behav_YA(:,1), brainscore_YA, 'rows', 'complete');
    if r_stable_YA > 0
        stat = flip_behav_pls_lv(stat, lv);
    end
end
end

function stat = flip_behav_pls_lv(stat, lv)
stat.brainscores{1}(:,lv) = -stat.brainscores{1}(:,lv);
stat.brainscores{2}(:,lv) = -stat.brainscores{2}(:,lv);
stat.behavscores{1}(:,lv) = -stat.behavscores{1}(:,lv);
stat.behavscores{2}(:,lv) = -stat.behavscores{2}(:,lv);

if isfield(stat, 'boot_res') && isfield(stat.boot_res, 'ulcorr') && isfield(stat.boot_res, 'llcorr')
    tmp_ul = stat.boot_res.ulcorr(:,lv);
    tmp_ll = stat.boot_res.llcorr(:,lv);
    stat.boot_res.ulcorr(:,lv) = -tmp_ll;
    stat.boot_res.llcorr(:,lv) = -tmp_ul;
end
end

function plot_lv_map(fig, panel_pos, stat_in, lv, topo_windows, cmap, panel_title)
stat = stat_in;
statdims = size(stat.stat);
bsr = reshape(stat.results.boot_result.compare_u(:,lv), statdims);
stat.stat = bsr;

if isempty(topo_windows)
    topo_windows(1).timerange = [min(stat.time) max(stat.time)];
    topo_windows(1).freqrange = [];
    topo_windows(1).label = 'All time';
end

if lv == 1
    clim_tfr = [-3 3];
    clim_topo = [-6 6];
else
    clim_tfr = [-1.5 1.5];
    clim_topo = [-3 3];
end

left = panel_pos(1);
bottom = panel_pos(2);
width = panel_pos(3);
height = panel_pos(4);

titleAx = axes(fig, 'Position', [left bottom + height * 0.910 width height * 0.080]);
axis(titleAx, 'off')
text(titleAx, 0, 0.5, panel_title, 'FontWeight', 'bold', 'FontSize', 9, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle');

nTopo = numel(topo_windows);
topoW = min(width * 0.280, width / max(nTopo, 2) * 0.55);
topoGap = width * 0.018;
topoTotalW = nTopo * topoW + (nTopo - 1) * topoGap;
topoStart = left + width * 0.16 + (width * 0.50 - topoTotalW) / 2;
if nTopo > 1
    topoStart = left + width * 0.205;
end
topoY = bottom + height * 0.545;
topoH = height * 0.300;

s = stat.results.s;
p_lv = stat.results.perm_result.sprob(lv);
r_lv = latent_lv_corr(stat, lv);
p_txt = sprintf('p = %.3f', p_lv);
if p_lv < .001
    p_txt = 'p < .001';
end
statsAxW = max(topoStart - left - width * 0.020, width * 0.120);
statsAx = axes(fig, 'Position', [left + width * 0.005, topoY + topoH * 0.10, statsAxW, topoH * 0.70]);
axis(statsAx, 'off')
text(statsAx, 0, 0.62, sprintf('LV%d\n%s\nExpl = %.2f\nr = %.2f', ...
    lv, p_txt, s(lv).^2 / sum(s.^2), r_lv), ...
    'Units', 'normalized', 'FontSize', 5.0, 'VerticalAlignment', 'middle', ...
    'HorizontalAlignment', 'left')

for iw = 1:nTopo
    topo = make_topo(stat, bsr, topo_windows(iw).timerange, topo_windows(iw).freqrange);
    topoAx = axes(fig, 'Position', [topoStart + (iw - 1) * (topoW + topoGap), topoY, topoW, topoH]);
    plot_topo(topoAx, topo, clim_topo, cmap);
    text(topoAx, 0.5, -0.11, topo_windows(iw).label, 'Units', 'normalized', ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', ...
        'FontSize', 5.5, 'FontWeight', 'normal')
end

cbTopoAx = axes(fig, 'Position', [topoStart + topoTotalW + width * 0.018, topoY + topoH * 0.18, width * 0.014, topoH * 0.58]);
imagesc(cbTopoAx, [0 1], linspace(clim_topo(1), clim_topo(2), 256)', linspace(clim_topo(1), clim_topo(2), 256)');
set(cbTopoAx, 'YDir', 'normal', 'XTick', [], 'YAxisLocation', 'right', ...
    'YTick', [clim_topo(1) 0 clim_topo(2)], 'FontSize', 6, 'TickDir', 'out')
ylabel(cbTopoAx, 'BSR', 'FontSize', 4.2)
colormap(cbTopoAx, cmap)
box(cbTopoAx, 'off')

tfrAx = axes(fig, 'Position', [left + width * 0.11, bottom + height * 0.115, width * 0.68, height * 0.330]);
tfr_bsr = make_tfr_map(stat, bsr);
imagesc(tfrAx, stat.time, stat.freq, tfr_bsr);
set(tfrAx, 'YDir', 'normal')
colormap(tfrAx, cmap)
caxis(tfrAx, clim_tfr)
hold(tfrAx, 'on')
xline(tfrAx, 0, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off')
for iw = 1:numel(topo_windows)
    tr = topo_windows(iw).timerange;
    fr = topo_windows(iw).freqrange;
    if isempty(fr)
        fr = [min(stat.freq) max(stat.freq)];
    end
    plot(tfrAx, [tr(1) tr(2) tr(2) tr(1) tr(1)], ...
        [fr(1) fr(1) fr(2) fr(2) fr(1)], ...
        'k-', 'LineWidth', 0.45)
end
xlabel(tfrAx, 'Time from blink peak (s)', 'FontSize', 7)
ylabel(tfrAx, 'Time scale (ms)', 'FontSize', 7)
box(tfrAx, 'off')
tfrAx.FontSize = 7;
tfrAx.TickDir = 'out';
tfrAx.XTick = [-1 0 1];
tfrAx.YTick = stat.freq;
dt = median(diff(stat.time));
df = median(diff(stat.freq));
tfrAx.XLim = [min(stat.time)-dt/2 max(stat.time)+dt/2];
tfrAx.YLim = [min(stat.freq)-df/2 max(stat.freq)+df/2];

cb = colorbar(tfrAx);
cb.Box = 'off';
cb.TickDirection = 'out';
cb.FontSize = 6;
cb.Ticks = [clim_tfr(1) 0 clim_tfr(2)];
cb.Label.String = 'BSR';
cb.Label.FontSize = 4.2;
cb.Position = [left + width * 0.82, bottom + height * 0.175, width * 0.014, height * 0.205];
end

function r = latent_lv_corr(stat, lv)
brain = [zscore(stat.brainscores{1}(:,lv)); zscore(stat.brainscores{2}(:,lv))];
behav = [zscore(stat.behavscores{1}(:,lv)); zscore(stat.behavscores{2}(:,lv))];
r = corr(brain, behav, 'rows', 'complete', 'Type', 'Pearson');
end

function topo = make_topo(stat, bsr, timerange, freqrange)
time_sel = stat.time >= timerange(1) & stat.time <= timerange(2);
if isempty(freqrange)
    freq_sel = true(size(stat.freq));
else
    freq_sel = stat.freq >= freqrange(1) & stat.freq <= freqrange(2);
end
topo_avg = squeeze(mean(mean(bsr(:, freq_sel, time_sel), 3, 'omitnan'), 2, 'omitnan'));
topo = [];
topo.label = stat.label;
topo.avg = topo_avg;
topo.dimord = 'chan';
topo.time = 0;
end

function tfr_bsr = make_tfr_map(stat, bsr)
tfr_bsr = squeeze(mean(bsr, 1, 'omitnan'));
end

function plot_topo(ax, topo, clim_topo, cmap)
axes(ax)
cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.parameter = 'avg';
cfg.colorbar = 'no';
cfg.comment = 'no';
cfg.marker = 'off';
cfg.figure = 'gca';
cfg.gridscale = 180;
cfg.style = 'both';
cfg.contournum = 3;
cfg.colormap = cmap;
ft_topoplotER(cfg, topo);
colormap(ax, cmap)
caxis(ax, clim_topo)
axis(ax, 'tight')
h = findall(ax, 'Type', 'line');
set(h, 'LineWidth', 0.24);
end

function plot_behav_corr_panel(fig, panel_pos, stat, lv, colYAbar, colOAbar, colYAerr, colOAerr, panel_title, show_legend)
[r_YA, r_OA, r_ALL, errY_low, errY_high, errO_low, errO_high, errA_low, errA_high, labels, sig_labels] = ...
    compute_domain_correlations(stat, lv);

ax = axes(fig, 'Position', panel_pos);
hold(ax, 'on')

R = [r_YA; r_OA]';
b = bar(ax, R, 'grouped', 'BarWidth', 0.75);
b(1).FaceColor = colYAbar;
b(1).EdgeColor = 'none';
b(2).FaceColor = colOAbar;
b(2).EdgeColor = 'none';

xYA = b(1).XEndPoints;
xOA = b(2).XEndPoints;

errorbar(ax, xYA, r_YA, errY_low, errY_high, ...
    'Color', colYAerr, 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 6);
errorbar(ax, xOA, r_OA, errO_low, errO_high, ...
    'Color', colOAerr, 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 6);

xMid = 1:numel(labels);
hAll = errorbar(ax, xMid, r_ALL, errA_low, errA_high, ...
    '-o', 'Color', [0 0 0], 'MarkerFaceColor', [0 0 0], ...
    'MarkerEdgeColor', [0 0 0], 'LineWidth', 1.0, ...
    'MarkerSize', 3.5, 'CapSize', 6);

for i = 1:numel(r_YA)
    if sig_labels(i) ~= ""
        ysig = max([r_YA(i)+errY_high(i), r_OA(i)+errO_high(i), r_ALL(i)+errA_high(i)]) + 0.05;
        text(ax, mean([xYA(i), xOA(i)]), ysig, sig_labels(i), ...
            'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold')
    end
end

yline(ax, 0, 'k-', 'LineWidth', 0.65)
ylim(ax, [-1 1])
xlim(ax, [0.5 numel(labels)+0.5])
xticks(ax, 1:numel(labels))
xticklabels(ax, labels)
xtickangle(ax, 32)
ylabel(ax, 'Pearson''s r')
title(ax, panel_title, 'FontWeight', 'bold')
if show_legend
    legend(ax, [b(1) b(2) hAll], {'Young','Older','Partial r | group'}, 'Location', 'southoutside', ...
        'Orientation', 'horizontal', 'Box', 'off')
end
format_axes(ax)
end

function [r_YA, r_OA, r_ALL, errY_low, errY_high, errO_low, errO_high, errA_low, errA_high, labels, sig_labels] = compute_domain_correlations(stat, lv)
corrtype = 'Pearson';
labels = {'Cryst.','Attention','Learning','Working','Fluid'};

% Match plot_PLSCcorr_blink.m: correlate brain scores with the exact
% behavioral design matrix used by the PLS, not a separately reconstructed
% and IQR-cleaned raw behavior table.
ord = [5 4 3 1 2];
behav_YA = stat.results.stacked_behavdata(1:size(stat.brainscores{1},1), ord);
behav_OA = stat.results.stacked_behavdata(size(stat.brainscores{1},1)+1:end, ord);

brainscore_YA = stat.brainscores{1}(:,lv);
brainscore_OA = stat.brainscores{2}(:,lv);
brainscore_ALL = [brainscore_YA; brainscore_OA];
behav_ALL = [behav_YA; behav_OA];
group_ALL = [zeros(size(brainscore_YA)); ones(size(brainscore_OA))];

r_YA = nan(1, numel(labels));
r_OA = nan(1, numel(labels));
r_ALL = nan(1, numel(labels));
for ib = 1:numel(labels)
    r_YA(ib) = corr(behav_YA(:,ib), brainscore_YA, 'rows', 'complete', 'Type', corrtype);
    r_OA(ib) = corr(behav_OA(:,ib), brainscore_OA, 'rows', 'complete', 'Type', corrtype);
    validA = ~isnan(behav_ALL(:,ib)) & ~isnan(brainscore_ALL) & ~isnan(group_ALL);
    r_ALL(ib) = partialcorr(behav_ALL(validA,ib), brainscore_ALL(validA), group_ALL(validA), 'Type', corrtype);
end

% Use the PLS result confidence intervals for the bar uncertainty. The
% plotted bar heights are raw-z manifest correlations, so transfer CI width
% rather than mixing in a second, noisier bootstrap over cleaned raw domains.
ll = stat.boot_res.llcorr;
ul = stat.boot_res.ulcorr;

ciw_YA = (ul(ord,lv)' - ll(ord,lv)') ./ 2;
ciw_OA = (ul(5+ord,lv)' - ll(5+ord,lv)') ./ 2;

errY_low = ciw_YA;
errY_high = ciw_YA;
errO_low = ciw_OA;
errO_high = ciw_OA;

rng(1)
nBoot = 5000;
ll_ALL = nan(size(r_ALL));
ul_ALL = nan(size(r_ALL));
for ib = 1:numel(labels)
    [ll_ALL(ib), ul_ALL(ib)] = bootstrap_partialcorr_ci(behav_ALL(:,ib), brainscore_ALL, group_ALL, corrtype, nBoot);
end
errA_low = r_ALL - ll_ALL;
errA_high = ul_ALL - r_ALL;

sig_labels = strings(1, numel(labels));
p_diff = nan(1, numel(labels));
for i = 1:numel(labels)
    validY = ~isnan(behav_YA(:,i)) & ~isnan(brainscore_YA);
    validO = ~isnan(behav_OA(:,i)) & ~isnan(brainscore_OA);
    nY = sum(validY);
    nO = sum(validO);
    zY = atanh(r_YA(i));
    zO = atanh(r_OA(i));
    se = sqrt(1/(nY-3) + 1/(nO-3));
    p_diff(i) = 2 * (1 - normcdf(abs((zO - zY) / se))); %#ok<AGROW>
end
p_fdr = mafdr(p_diff, 'BHFDR', true);
for i = 1:numel(labels)
    if p_fdr(i) < .001
        sig_labels(i) = "***";
    elseif p_fdr(i) < .01
        sig_labels(i) = "**";
    elseif p_fdr(i) < .05
        sig_labels(i) = "*";
    end
end
end

function [ll, ul] = bootstrap_partialcorr_ci(x, y, g, corrtype, nBoot)
valid = ~isnan(x) & ~isnan(y) & ~isnan(g);
x = x(valid);
y = y(valid);
g = g(valid);
n = numel(x);
boot_r = nan(nBoot,1);
for iboot = 1:nBoot
    idx = randi(n, n, 1);
    boot_r(iboot) = partialcorr(x(idx), y(idx), g(idx), 'Type', corrtype);
end
ll = prctile(boot_r, 2.5);
ul = prctile(boot_r, 97.5);
end

function plot_lv1_lv3_scatter_panel(fig, panel_pos, stat, colYA, colOA)
lv1_YA = stat.brainscores{1}(:,1);
lv3_YA = stat.brainscores{1}(:,3);
lv1_OA = stat.brainscores{2}(:,1);
lv3_OA = stat.brainscores{2}(:,3);

[r_YA, p_YA] = corr(lv1_YA, lv3_YA, 'rows', 'complete', 'Type', 'Pearson');
[r_OA, p_OA] = corr(lv1_OA, lv3_OA, 'rows', 'complete', 'Type', 'Pearson');
[r_ALL, p_ALL] = corr([lv1_YA; lv1_OA], [lv3_YA; lv3_OA], ...
    'rows', 'complete', 'Type', 'Pearson');

lv1_all_z = zscore([lv1_YA; lv1_OA]);
lv3_all_z = zscore([lv3_YA; lv3_OA]);
lv1_YA = lv1_all_z(1:numel(lv1_YA));
lv1_OA = lv1_all_z(numel(lv1_YA)+1:end);
lv3_YA = lv3_all_z(1:numel(lv3_YA));
lv3_OA = lv3_all_z(numel(lv3_YA)+1:end);

ax = axes(fig, 'Position', panel_pos);
hold(ax, 'on')
scatter(ax, lv1_YA, lv3_YA, 14, colYA, 'filled', ...
    'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'w', 'LineWidth', 0.35);
scatter(ax, lv1_OA, lv3_OA, 14, colOA, 'filled', ...
    'MarkerFaceAlpha', 0.65, 'MarkerEdgeColor', 'w', 'LineWidth', 0.35);
plot_lm_line(ax, lv1_YA, lv3_YA, colYA, '-', 1.0);
plot_lm_line(ax, lv1_OA, lv3_OA, colOA, '-', 1.0);
plot_lm_line(ax, [lv1_YA; lv1_OA], [lv3_YA; lv3_OA], [0.10 0.10 0.10], ':', 1.2);

xlabel(ax, 'LV1 brain score')
ylabel(ax, 'LV3 brain score')
title(ax, 'E  LV1 vs LV3 brain scores', 'FontWeight', 'bold')
txt = sprintf('Y r=%.2f, %s\nO r=%.2f, %s\nAll r=%.2f, %s', ...
    r_YA, p_string(p_YA), r_OA, p_string(p_OA), r_ALL, p_string(p_ALL));
text(ax, 0.03, 0.97, txt, 'Units', 'normalized', 'FontSize', 4.1, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
axis(ax, 'padded')
format_axes(ax)
end

function plot_lm_line(ax, x, y, col, style, line_width)
valid = ~isnan(x) & ~isnan(y);
if sum(valid) <= 2
    return
end
x = x(valid);
y = y(valid);
pfit = polyfit(x, y, 1);
xfit = linspace(min(x), max(x), 100);
plot(ax, xfit, polyval(pfit, xfit), 'Color', col, ...
    'LineStyle', style, 'LineWidth', line_width)
end

function plot_interaction_panel(fig, panel_pos, results, colYA, colOA)
regime_names = {'Stable','Balanced','Flexible'};
idx_regime = results.Type == "Regime";
regime_results = results(idx_regime,:);
[~, plot_ord] = ismember(regime_names, cellstr(regime_results.Outcome));
if any(plot_ord == 0)
    error('Could not find all regime rows in LV1 x LV3 results table.');
end

x = 1:3;
yYA = regime_results.B_LV1xLV3_YA(plot_ord);
yOA = regime_results.B_LV1xLV3_OA(plot_ord);
eYA = regime_results.SE_LV1xLV3_YA(plot_ord);
eOA = regime_results.SE_LV1xLV3_OA(plot_ord);
pYA_trend = endpoint_trend_p(yYA, eYA);
pOA_trend = endpoint_trend_p(yOA, eOA);
p_group = endpoint_group_trend_p(yYA, eYA, yOA, eOA);

ax = axes(fig, 'Position', panel_pos);
hold(ax, 'on')

fill(ax, [x fliplr(x)], [yYA' + eYA' fliplr(yYA' - eYA')], colYA, ...
    'FaceAlpha', 0.16, 'EdgeColor', 'none');
fill(ax, [x fliplr(x)], [yOA' + eOA' fliplr(yOA' - eOA')], colOA, ...
    'FaceAlpha', 0.16, 'EdgeColor', 'none');
hYA = plot(ax, x, yYA, 'o-', 'Color', colYA, 'MarkerFaceColor', colYA, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.6, 'MarkerSize', 5);
hOA = plot(ax, x, yOA, 'o-', 'Color', colOA, 'MarkerFaceColor', colOA, ...
    'MarkerEdgeColor', 'w', 'LineWidth', 1.6, 'MarkerSize', 5);
errorbar(ax, x, yYA, eYA, 'Color', colYA, 'LineStyle', 'none', 'LineWidth', 0.9, 'CapSize', 6);
errorbar(ax, x, yOA, eOA, 'Color', colOA, 'LineStyle', 'none', 'LineWidth', 0.9, 'CapSize', 6);

yline(ax, 0, 'k-', 'LineWidth', 0.75)
xlim(ax, [0.75 3.36])
xticks(ax, x)
xticklabels(ax, regime_names)
ylabel(ax, 'LV1 x LV3 coefficient')
xlabel(ax, 'Cognitive regime')
yl = ylim(ax);
title(ax, 'F  LV1 x LV3 crossover', 'FontWeight', 'bold')
text(ax, 3.04, yYA(end), 'Young', 'FontSize', 6, 'Color', colYA, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
text(ax, 3.04, yOA(end), 'Older', 'FontSize', 6, 'Color', colOA, ...
    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle')
txt = sprintf('Stable-Flexible trend:\nY %s\nO %s\nY vs O %s', ...
    p_string(pYA_trend), p_string(pOA_trend), p_string(p_group));
text(ax, 0.03, 0.97, txt, 'Units', 'normalized', 'FontSize', 4.1, ...
    'VerticalAlignment', 'top', 'HorizontalAlignment', 'left')
ylim(ax, yl)
format_axes(ax)
end

function p = endpoint_trend_p(y, se)
contrast = y(end) - y(1);
contrast_se = sqrt(se(end).^2 + se(1).^2);
z = contrast ./ contrast_se;
p = 2 .* (1 - normcdf(abs(z)));
end

function p = endpoint_group_trend_p(yYA, eYA, yOA, eOA)
contrast = (yYA(end) - yYA(1)) - (yOA(end) - yOA(1));
contrast_se = sqrt(eYA(end).^2 + eYA(1).^2 + eOA(end).^2 + eOA(1).^2);
z = contrast ./ contrast_se;
p = 2 .* (1 - normcdf(abs(z)));
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
set(ax, 'FontSize', 8, 'LineWidth', 0.75, 'TickDir', 'out', ...
    'XColor', 'k', 'YColor', 'k');
end
