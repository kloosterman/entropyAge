%% Correlate LV1 and LV3 brain scores
% Scatter plot with pooled and group-specific regression lines.

folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

if ~exist(plotfolder, 'dir')
    mkdir(plotfolder);
end

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
                fprintf('Loaded stat_mse_2group_behav from %s\n', char(candidate_files(ifile)));
                break
            elseif isfield(S, 'stat_mse_2group')
                stat_mse_2group_behav = S.stat_mse_2group;
                loaded_stat = true;
                fprintf('Loaded stat_mse_2group as stat_mse_2group_behav from %s\n', char(candidate_files(ifile)));
                break
            end
        end
    end

    if ~loaded_stat
        error(['stat_mse_2group_behav is not in the workspace and no fallback file was found. ', ...
            'Run runPLSanalyses first or save the stat struct in stats_structs.']);
    end
end

lv_x = 1;
lv_y = 3;
corrtype = 'Pearson';

brainscore_YA_LV1 = stat_mse_2group_behav.brainscores{1}(:,lv_x);
brainscore_YA_LV3 = stat_mse_2group_behav.brainscores{1}(:,lv_y);

brainscore_OA_LV1 = stat_mse_2group_behav.brainscores{2}(:,lv_x);
brainscore_OA_LV3 = stat_mse_2group_behav.brainscores{2}(:,lv_y);

brainscore_ALL_LV1 = [brainscore_YA_LV1; brainscore_OA_LV1];
brainscore_ALL_LV3 = [brainscore_YA_LV3; brainscore_OA_LV3];
group_ALL = [repmat("YA", size(brainscore_YA_LV1)); repmat("OA", size(brainscore_OA_LV1))];

[r_YA, p_YA] = corr(brainscore_YA_LV1, brainscore_YA_LV3, ...
    'rows', 'complete', 'Type', corrtype);
[r_OA, p_OA] = corr(brainscore_OA_LV1, brainscore_OA_LV3, ...
    'rows', 'complete', 'Type', corrtype);
[r_ALL, p_ALL] = corr(brainscore_ALL_LV1, brainscore_ALL_LV3, ...
    'rows', 'complete', 'Type', corrtype);

n_YA = sum(~isnan(brainscore_YA_LV1) & ~isnan(brainscore_YA_LV3));
n_OA = sum(~isnan(brainscore_OA_LV1) & ~isnan(brainscore_OA_LV3));
n_ALL = sum(~isnan(brainscore_ALL_LV1) & ~isnan(brainscore_ALL_LV3));

fprintf('\nLV%d vs LV%d brain score correlations (%s)\n', lv_x, lv_y, corrtype);
fprintf('YA:  r = %.3f, p = %.4f, n = %d\n', r_YA, p_YA, n_YA);
fprintf('OA:  r = %.3f, p = %.4f, n = %d\n', r_OA, p_OA, n_OA);
fprintf('All: r = %.3f, p = %.4f, n = %d\n', r_ALL, p_ALL, n_ALL);

results_table = table( ...
    ["YA"; "OA"; "All"], ...
    [r_YA; r_OA; r_ALL], ...
    [p_YA; p_OA; p_ALL], ...
    [n_YA; n_OA; n_ALL], ...
    'VariableNames', {'Group', 'r', 'p', 'n'});

col_YA = [0.85 0.25 0.25];
col_OA = [0.20 0.45 0.85];
col_ALL = [0.15 0.15 0.15];

f = figure('Color', 'w');
f.Position = [200 200 520 460];
hold on

hYA = scatter(brainscore_YA_LV1, brainscore_YA_LV3, ...
    38, col_YA, 'filled', ...
    'MarkerFaceAlpha', 0.70, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.5);

hOA = scatter(brainscore_OA_LV1, brainscore_OA_LV3, ...
    38, col_OA, 'filled', ...
    'MarkerFaceAlpha', 0.70, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.5);

valid_YA = ~isnan(brainscore_YA_LV1) & ~isnan(brainscore_YA_LV3);
if sum(valid_YA) > 2
    x = brainscore_YA_LV1(valid_YA);
    y = brainscore_YA_LV3(valid_YA);
    pfit = polyfit(x, y, 1);
    xfit = linspace(min(x), max(x), 200);
    hFitYA = plot(xfit, polyval(pfit, xfit), ...
        'Color', col_YA, ...
        'LineWidth', 1.6);
else
    hFitYA = gobjects(1);
end

valid_OA = ~isnan(brainscore_OA_LV1) & ~isnan(brainscore_OA_LV3);
if sum(valid_OA) > 2
    x = brainscore_OA_LV1(valid_OA);
    y = brainscore_OA_LV3(valid_OA);
    pfit = polyfit(x, y, 1);
    xfit = linspace(min(x), max(x), 200);
    hFitOA = plot(xfit, polyval(pfit, xfit), ...
        'Color', col_OA, ...
        'LineWidth', 1.6);
else
    hFitOA = gobjects(1);
end

valid_ALL = ~isnan(brainscore_ALL_LV1) & ~isnan(brainscore_ALL_LV3);
if sum(valid_ALL) > 2
    x = brainscore_ALL_LV1(valid_ALL);
    y = brainscore_ALL_LV3(valid_ALL);
    pfit = polyfit(x, y, 1);
    xfit = linspace(min(x), max(x), 200);
    hFitALL = plot(xfit, polyval(pfit, xfit), ...
        'Color', col_ALL, ...
        'LineStyle', '--', ...
        'LineWidth', 2.0);
else
    hFitALL = gobjects(1);
end

xlabel(sprintf('LV%d brain score', lv_x));
ylabel(sprintf('LV%d brain score', lv_y));
title(sprintf('LV%d vs LV%d brain scores', lv_x, lv_y));

legend([hYA hOA hFitYA hFitOA hFitALL], ...
    {sprintf('YA: r=%.2f, p=%.3f', r_YA, p_YA), ...
     sprintf('OA: r=%.2f, p=%.3f', r_OA, p_OA), ...
     'YA fit', ...
     'OA fit', ...
     sprintf('All fit: r=%.2f, p=%.3f', r_ALL, p_ALL)}, ...
    'Location', 'best', ...
    'Box', 'off');

axis padded
box on
grid on
set(gca, ...
    'FontSize', 9, ...
    'LineWidth', 1, ...
    'TickDir', 'out', ...
    'GridAlpha', 0.12);

outfile = sprintf('LV%d_vs_LV%d_brainscore_scatter', lv_x, lv_y);

if exist('exportFigure', 'file')
    exportFigure(f, fullfile(plotfolder, [outfile '.pdf']), 'FontSize', 9);
    exportFigure(f, fullfile(plotfolder, [outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f, fullfile(plotfolder, [outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f, fullfile(plotfolder, [outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

writetable(results_table, fullfile(plotfolder, [outfile '_corrs.csv']));
save(fullfile(plotfolder, [outfile '_corrs.mat']), ...
    'results_table', ...
    'brainscore_YA_LV1', 'brainscore_YA_LV3', ...
    'brainscore_OA_LV1', 'brainscore_OA_LV3', ...
    'brainscore_ALL_LV1', 'brainscore_ALL_LV3', ...
    'group_ALL');

fprintf('Saved %s outputs to %s\n', outfile, char(plotfolder));
