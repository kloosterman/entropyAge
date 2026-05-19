%% plot correlations task PLS brainscores with raw behavior
% Figure
plotinfig = 1
if ~plotinfig
    f = figure;

    % auto layout depending on number of panels
    plotWhat = 'domains';   % 'regimes' or 'domains'

    if strcmpi(plotWhat,'regimes')
        tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
    else
        tiledlayout(3,4,'TileSpacing','compact','Padding','compact');
    end
end

figLayout = findobj(gcf, 'Type', 'tiledlayout');
if isempty(figLayout)
    figLayout = tiledlayout(3,6,'TileSpacing','compact','Padding','tight');
else
    figLayout = figLayout(1);
end

% =========================
% SETTINGS
% =========================
useBinnedAll = 0;        % plot binned across-subject average
nBins        = 5;
showRawPts   = true;
showBinErr   = true;
minPerBin    = 3;
corrtype     = 'Pearson';

% -----------------------------------------
% SUBJECT EXCLUSION: define once
% indices are within group
% -----------------------------------------
ya_excl = [];    % e.g. [3 8]
oa_excl = 5;     % e.g. [7]

% Data
Y_BStask = stat_mse_2group_task.brainscores{1}(:,1);
O_BStask = stat_mse_2group_task.brainscores{2}(:,1);

nYA_full = length(Y_BStask);
nOA_full = length(O_BStask);

% =========================
% CHOOSE WHAT TO PLOT
% =========================
clear lower
switch lower(plotWhat)
    case 'regimes'
        behavdat_full = table2array(regime_tbl);
        plot_labels = regime_tbl.Properties.VariableNames;
        behavoi = 1:size(behavdat_full,2);
        behavoi = fliplr(behavoi); % start with stable

        % tile positions for 3-panel layout
        tile_idx = 7:2:11;

    case 'domains'
        domainNames = {'Crystallized','Attention','LearningMemory','WorkingMemory','FluidIntelligence'};
        plot_labels = {'Crystallized','Attention','Learning & Memory','Working Memory','Fluid Intelligence'};
        behavdat_full = table2array(clean_tbl(:, domainNames));
        behavoi = 1:numel(domainNames);

        % two-column panels in rows 2-3 of the 3x6 master layout
        tile_idx = [7 9 11 13 15];
end

% =========================
% APPLY SAME EXCLUSION TO BRAIN + BEHAVIOR
% =========================
keepYA = true(nYA_full,1);
keepOA = true(nOA_full,1);

keepYA(ya_excl) = false;
keepOA(oa_excl) = false;

Y_BStask = Y_BStask(keepYA);
O_BStask = O_BStask(keepOA);

behav_YA = behavdat_full(1:nYA_full, :);
behav_OA = behavdat_full(nYA_full+1:end, :);

behav_YA = behav_YA(keepYA, :);
behav_OA = behav_OA(keepOA, :);

% recombine after matching exclusions
BStask   = [Y_BStask; O_BStask];
behavdat = [behav_YA; behav_OA];

% z-score brain scores after exclusion
BStask = zscore(BStask);

nYA = sum(keepYA);
nOA = sum(keepOA);

Y_BStask = BStask(1:nYA);
O_BStask = BStask(nYA+1:end);

% Common axis limits
x_all = BStask;
xlim_all = [min(x_all) max(x_all)];
ylim_all = [min(behavdat(:)) max(behavdat(:))];

cols = [0.85 0.20 0.20; ...
        0.10 0.45 0.75];

for ii = 1:numel(behavoi)
    i = behavoi(ii);

    nexttile(figLayout, tile_idx(ii), [1 2]);
    hold on

    % -------------------------
    % Scatter
    % -------------------------
    if showRawPts
        h1 = scatter(Y_BStask, behav_YA(:,i), 15, 'filled', ...
            'MarkerFaceColor', cols(1,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.25);

        h2 = scatter(O_BStask, behav_OA(:,i), 15, 'filled', ...
            'MarkerFaceColor', cols(2,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.25);
    end

    % -------------------------
    % Linear fits within group
    % -------------------------
    pY_lin = polyfit(Y_BStask, behav_YA(:,i), 1);
    pO_lin = polyfit(O_BStask, behav_OA(:,i), 1);

    xY = linspace(min(Y_BStask), max(Y_BStask), 100);
    xO = linspace(min(O_BStask), max(O_BStask), 100);

    plot(xY, polyval(pY_lin, xY), '-', 'Color', cols(1,:), 'LineWidth', 1)
    plot(xO, polyval(pO_lin, xO), '-', 'Color', cols(2,:), 'LineWidth', 1)

    % -------------------------
    % Correlations
    % -------------------------
    [rY, pY] = corr(Y_BStask, behav_YA(:,i), 'rows', 'complete', 'Type', corrtype);
    [rO, pO] = corr(O_BStask, behav_OA(:,i), 'rows', 'complete', 'Type', corrtype);

    y_all_i = behavdat(:,i);
    [r_all, p_all] = corr(BStask, y_all_i, 'rows', 'complete', 'Type', corrtype);

    % -------------------------
    % Across-subject trend + model comparison
    % -------------------------
    xFit = linspace(min(BStask), max(BStask), 200);

    mdl_lin  = fitlm(BStask, y_all_i);
    mdl_quad = fitlm(BStask, y_all_i, 'quadratic');

    SSE_lin  = mdl_lin.SSE;
    SSE_quad = mdl_quad.SSE;

    df_lin   = mdl_lin.DFE;
    df_quad  = mdl_quad.DFE;

    df1      = df_lin - df_quad;
    df2      = df_quad;

    F_quad   = ((SSE_lin - SSE_quad) / df1) / (SSE_quad / df2);
    F_quad   = max(F_quad, 0);
    p_quad   = 1 - fcdf(F_quad, df1, df2);

    R2_lin   = mdl_lin.Rsquared.Adjusted;
    R2_quad  = mdl_quad.Rsquared.Adjusted;
    dR2      = R2_quad - R2_lin;

    if useBinnedAll
        edges = linspace(min(BStask), max(BStask), nBins+1);

        yBinMean = nan(1,nBins);
        yBinSEM  = nan(1,nBins);
        xBinMean = nan(1,nBins);
        nPerBin  = nan(1,nBins);

        binIdx = discretize(BStask, edges);
        binIdx(BStask == edges(end)) = nBins;

        for b = 1:nBins
            idxb = binIdx == b;
            nPerBin(b) = sum(idxb);

            if nPerBin(b) >= minPerBin
                xBinMean(b) = mean(BStask(idxb), 'omitnan');
                ytmp = y_all_i(idxb);
                yBinMean(b) = mean(ytmp, 'omitnan');
                yBinSEM(b)  = std(ytmp, 'omitnan') / sqrt(sum(~isnan(ytmp)));
            end
        end

        valid = ~isnan(xBinMean) & ~isnan(yBinMean);

        if showBinErr
            errorbar(xBinMean(valid), yBinMean(valid), yBinSEM(valid), ...
                'o-', 'Color', 'k', 'LineWidth', 1.8, ...
                'MarkerFaceColor', 'k', 'MarkerSize', 5, ...
                'CapSize', 0);
        else
            plot(xBinMean(valid), yBinMean(valid), 'o-', ...
                'Color', 'k', 'LineWidth', 1.8, ...
                'MarkerFaceColor', 'k', 'MarkerSize', 5);
        end

        if sum(valid) >= 4
            if strcmpi(plotWhat,'regimes') && i == 2
                pBin = polyfit(xBinMean(valid), yBinMean(valid), 2);
            else
                pBin = polyfit(xBinMean(valid), yBinMean(valid), 1);
            end
            yFitBin = polyval(pBin, xFit);
            plot(xFit, yFitBin, 'k--', 'LineWidth', 1.4)
        end

    else
        % if strcmpi(plotWhat,'regimes') && i == 2
        %     yFit = predict(mdl_quad, xFit');
        % else
            yFit = predict(mdl_lin, xFit');
        % end
        plot(xFit, yFit, 'Color', [0.4 0.4 0.4], 'LineWidth', 1) 
    end

    % -------------------------
    % Print stats to command window
    % -------------------------
    fprintf('\n=== %s ===\n', plot_labels{i});
    fprintf('r_YA = %.3f, p = %.4f\n', rY, pY);
    fprintf('r_OA = %.3f, p = %.4f\n', rO, pO);
    fprintf('r_all = %.3f, p = %.4f\n', r_all, p_all);
    fprintf('Linear adj R2 = %.3f\n', R2_lin);
    fprintf('Quadratic adj R2 = %.3f\n', R2_quad);
    fprintf('Delta R2 = %.3f\n', dR2);
    fprintf('Quadratic vs linear: F(%d,%d) = %.3f, p = %.4f\n', ...
        df1, df2, F_quad, p_quad);

    % -------------------------
    % Subplot title + stat line
    % -------------------------
    % format numbers without leading zero
    r_num = regexprep(sprintf('%.2f', r_all), '^(-?)0\.', '$1.');
    p_num = regexprep(sprintf('%.3f', p_all), '^(-?)0\.', '$1.');
    p_quad_num = regexprep(sprintf('%.3f', p_quad), '^(-?)0\.', '$1.');

    if strcmpi(plotWhat,'regimes')
        if i == 3   % stable
            if p_all < 0.001
                stat_str = sprintf('r = %s, p < .001', r_num);
            else
                stat_str = sprintf('r = %s, p = %s', r_num, p_num);
            end

        % elseif i == 2   % balanced
        %     if p_quad < 0.001
        %         stat_str = 'quad > lin, p < .001';
        %     else
        %         stat_str = sprintf('quad > lin, p = %s', p_quad_num);
        %     end

        % elseif i == 1   % flexible
        elseif i < 3   % flexible and balanced
            if p_all < 0.001
                stat_str = sprintf('r = %s, p < .001', r_num);
            else
                stat_str = sprintf('r = %s, p = %s', r_num, p_num);
            end

        else
            stat_str = '';
        end
    else
        if p_all < 0.001
            stat_str = sprintf('r = %s, p < .001', r_num);
        else
            stat_str = sprintf('r = %s, p = %s', r_num, p_num);
        end
    end

    title({['\it' plot_labels{i}], ['\rm' stat_str]}, ...
        'Interpreter', 'tex', ...
        'FontWeight', 'normal', ...
        'FontSize', 7)

    xlabel('Brain score (z)')
    ylabel('Behavior (z)')

    xlim(xlim_all)
    ylim(ylim_all)
    box off
    axis square
    axis padded
    ax = gca;
    ax.TickDir = 'out';
end

h = annotation('textbox', [0.0311    0.6250    1.0000    0.0500], ...
    'String', 'Brain score vs. behavior across cognitive regimes', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center', ...
    'FontWeight', 'bold', ...
    'FontSize', 8, ...
    'Color', 'k');
% plotedit on
