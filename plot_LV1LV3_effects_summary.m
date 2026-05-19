%% Summary figure for LV1/LV3 coordination effects
% Panels:
%   A) LV1 vs LV3 brain-score scatter
%   B) LV1 x LV3 interaction gradient across regimes
%   C) LV1 main effect across regimes
%   D) LV3 main effect across regimes

folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

if ~exist(plotfolder, 'dir')
    mkdir(plotfolder);
end

%% Load or compute results table

results_file = fullfile(plotfolder, 'LV1xLV3_predict_behavior_results.mat');

if exist(results_file, 'file')
    Sres = load(results_file, 'results');
    results = Sres.results;
else
    test_LV1xLV3_predict_behavior
end

needed_cols = {'B_LV1xLV3_YA','SE_LV1xLV3_YA','B_LV1_YA','SE_LV1_YA','B_LV3_YA','SE_LV3_YA'};
if ~all(ismember(needed_cols, results.Properties.VariableNames))
    fprintf('Saved results table is missing main-effect columns; rerunning analysis script.\n');
    test_LV1xLV3_predict_behavior
end

%% Load PLS stat struct for scatter panel

if ~exist('stat_mse_2group_behav', 'var')
    candidate_files = [
        fullfile(folder, 'behavPLS_2group_blink.mat')
        fullfile(folder, 'stat_mse_2group_blink.mat')
        fullfile(folder, 'stat_mse_2group_blink_modulation.mat')
    ];

    loaded_stat = false;
    for ifile = 1:numel(candidate_files)
        if exist(candidate_files(ifile), 'file')
            S = load(candidate_files(ifile));
            if isfield(S, 'stat_mse_2group_behav')
                stat_mse_2group_behav = S.stat_mse_2group_behav;
                loaded_stat = true;
                break
            elseif isfield(S, 'stat_mse_2group')
                stat_mse_2group_behav = S.stat_mse_2group;
                loaded_stat = true;
                break
            end
        end
    end

    if ~loaded_stat
        error(['stat_mse_2group_behav is not in the workspace and no fallback file was found. ', ...
            'Run runPLSanalyses first or save the stat struct in stats_structs.']);
    end
end

%% Shared plotting settings

col_YA = [0.85 0.25 0.25];
col_OA = [0.20 0.45 0.85];
col_ALL = [0.10 0.10 0.10];

regime_values = [-1 0 1];
regime_names = {'Stable','Balanced','Flexible'};
idx_regime = results.Type == "Regime";
regime_results = results(idx_regime,:);
[~, plot_ord] = ismember(regime_names, cellstr(regime_results.Outcome));
if any(plot_ord == 0)
    error('Could not find all regime rows in results table.');
end

%% Figure

f = figure('Color', 'w');
f.Units = 'centimeters';
f.Position = [2 2 18 15];

t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

%% A) LV1 vs LV3 scatter

axA = nexttile(t, 1);
hold(axA, 'on')

lv1_YA = stat_mse_2group_behav.brainscores{1}(:,1);
lv3_YA = stat_mse_2group_behav.brainscores{1}(:,3);
lv1_OA = stat_mse_2group_behav.brainscores{2}(:,1);
lv3_OA = stat_mse_2group_behav.brainscores{2}(:,3);

[r_YA, p_YA] = corr(lv1_YA, lv3_YA, 'rows', 'complete', 'Type', 'Pearson');
[r_OA, p_OA] = corr(lv1_OA, lv3_OA, 'rows', 'complete', 'Type', 'Pearson');

hYA = scatter(axA, lv1_YA, lv3_YA, 24, col_YA, 'filled', ...
    'MarkerFaceAlpha', 0.70, 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);
hOA = scatter(axA, lv1_OA, lv3_OA, 24, col_OA, 'filled', ...
    'MarkerFaceAlpha', 0.70, 'MarkerEdgeColor', 'w', 'LineWidth', 0.4);

plot_lm_line(axA, lv1_YA, lv3_YA, col_YA, '-', 1.4);
plot_lm_line(axA, lv1_OA, lv3_OA, col_OA, '-', 1.4);
plot_lm_line(axA, [lv1_YA; lv1_OA], [lv3_YA; lv3_OA], col_ALL, '--', 1.6);

xlabel(axA, 'LV1 brain score')
ylabel(axA, 'LV3 brain score')
title(axA, sprintf('A  LV1 vs LV3\nYA r=%.2f, OA r=%.2f', r_YA, r_OA), 'FontWeight', 'bold')
legend(axA, [hYA hOA], {'Young','Older'}, 'Location', 'best', 'Box', 'off')
format_axes(axA);

%% B) Interaction gradient

axB = nexttile(t, 2);
plot_effect_gradient(axB, regime_values, regime_names, regime_results, plot_ord, ...
    'B_LV1xLV3_YA', 'SE_LV1xLV3_YA', ...
    'B_LV1xLV3_OA', 'SE_LV1xLV3_OA', ...
    col_YA, col_OA, ...
    'LV1 x LV3 coefficient', ...
    'B  Interaction gradient');

%% C) LV1 main effect

axC = nexttile(t, 3);
plot_effect_gradient(axC, regime_values, regime_names, regime_results, plot_ord, ...
    'B_LV1_YA', 'SE_LV1_YA', ...
    'B_LV1_OA', 'SE_LV1_OA', ...
    col_YA, col_OA, ...
    'LV1 coefficient', ...
    'C  LV1 main effect');

%% D) LV3 main effect

axD = nexttile(t, 4);
plot_effect_gradient(axD, regime_values, regime_names, regime_results, plot_ord, ...
    'B_LV3_YA', 'SE_LV3_YA', ...
    'B_LV3_OA', 'SE_LV3_OA', ...
    col_YA, col_OA, ...
    'LV3 coefficient', ...
    'D  LV3 main effect');

linkaxes([axB axC axD], 'x')

outfile = 'LV1_LV3_effects_summary';

if exist('exportFigure', 'file')
    exportFigure(f, fullfile(plotfolder, [outfile '.pdf']), 'FontSize', 9);
    exportFigure(f, fullfile(plotfolder, [outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f, fullfile(plotfolder, [outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f, fullfile(plotfolder, [outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

fprintf('Saved %s outputs to %s\n', outfile, char(plotfolder));

%% Local helpers

function plot_lm_line(ax, x, y, col, style, width)
valid = ~isnan(x) & ~isnan(y);
if sum(valid) <= 2
    return
end
x = x(valid);
y = y(valid);
pfit = polyfit(x, y, 1);
xfit = linspace(min(x), max(x), 100);
plot(ax, xfit, polyval(pfit, xfit), ...
    'Color', col, 'LineStyle', style, 'LineWidth', width);
end

function plot_effect_gradient(ax, xvals, xlabels, tbl, ord, bYA, seYA, bOA, seOA, colYA, colOA, ylab, ttl)
hold(ax, 'on')

x = xvals(:);
yYA = tbl.(bYA)(ord);
yOA = tbl.(bOA)(ord);
eYA = tbl.(seYA)(ord);
eOA = tbl.(seOA)(ord);
bar_width = 0.26;
bar_offset = 0.14;

hYA = bar(ax, x - bar_offset, yYA, bar_width, ...
    'FaceColor', colYA, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.72);
hOA = bar(ax, x + bar_offset, yOA, bar_width, ...
    'FaceColor', colOA, ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.72);

errorbar(ax, x - bar_offset, yYA, eYA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);
errorbar(ax, x + bar_offset, yOA, eOA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);

plot(ax, x - bar_offset, yYA, 'o-', ...
    'Color', colYA, ...
    'MarkerFaceColor', colYA, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.3, ...
    'MarkerSize', 4.5);
plot(ax, x + bar_offset, yOA, 'o-', ...
    'Color', colOA, ...
    'MarkerFaceColor', colOA, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.3, ...
    'MarkerSize', 4.5);

yline(ax, 0, 'k-', 'LineWidth', 0.75);
xlim(ax, [-1.2 1.2])
xticks(ax, xvals)
xticklabels(ax, xlabels)
ylabel(ax, ylab)
title(ax, ttl, 'FontWeight', 'bold')
legend(ax, [hYA hOA], {'Young','Older'}, 'Location', 'best', 'Box', 'off')
format_axes(ax);
end

function format_axes(ax)
box(ax, 'on')
set(ax, ...
    'FontSize', 9, ...
    'LineWidth', 1, ...
    'TickDir', 'out', ...
    'XColor', 'k', ...
    'YColor', 'k');
end
