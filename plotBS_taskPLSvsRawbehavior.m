%% plot correlations task PLS brainscores with raw behavior
% Figure
if not(plotinfig)
  f = figure;
  tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
end

% =========================
% SETTINGS
% =========================
useBinnedAll = 0;        % plot binned across-subject average
nBins        = 5;        % e.g. 5 or 7 bins
showRawPts   = true;     % keep YA/OA scatter
showBinErr   = true;     % error bars for binned means
minPerBin    = 3;        % require at least this many subjects per bin

% Data
Y_BStask = stat_mse_2group_task.brainscores{1}(:,1);
O_BStask = stat_mse_2group_task.brainscores{2}(:,1);
BStask   = [Y_BStask; O_BStask];

behavoi = [1 2 3];
regime_labels = regime_tbl.Properties.VariableNames;

behavdat = table2array(regime_tbl);

nYA = length(Y_BStask);

% Identify OA rows in full dataset
idx_OA = (1:length(BStask)) > nYA;

% Work on full behavior
y_all = behavdat(:,3);

% Compute z within OA only
z_OA = zscore(y_all(idx_OA));

% Logical index of subjects to keep
keep_OA = true(sum(idx_OA),1);
keep_OA(z_OA < -3) = false;

% Build full keep index
keep = true(size(BStask));
keep(idx_OA) = keep_OA;

% Apply to everything
BStask = BStask(keep);
behavdat = behavdat(keep,:);

% Re-split YA / OA after removal
nYA_new = sum(~idx_OA(keep));

BStask = zscore(BStask);

Y_BStask = BStask(1:nYA_new);
O_BStask = BStask(nYA_new+1:end);

behav_YA = behavdat(1:nYA_new,:);
behav_OA = behavdat(nYA_new+1:end,:);

% Common axis limits
x_all = BStask;
xlim_all = [min(x_all) max(x_all)];
ylim_all = [min(behavdat(:)) max(behavdat(:))];

cols = [1 0 0;          % Young = red
        0 0.447 0.741]; % Older = blue

for i = 1:3
    nexttile(i+3);
    hold on

    % -------------------------
    % Scatter
    % -------------------------
    if showRawPts
        h1 = scatter(Y_BStask, behav_YA(:,i), 30, 'filled', ...
            'MarkerFaceColor', cols(1,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.5);

        h2 = scatter(O_BStask, behav_OA(:,i), 30, 'filled', ...
            'MarkerFaceColor', cols(2,:), ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 0.5);
    end

    % -------------------------
    % Linear fits within group
    % -------------------------
    pY_lin = polyfit(Y_BStask, behav_YA(:,i), 1);
    pO_lin = polyfit(O_BStask, behav_OA(:,i), 1);

    xY = linspace(min(Y_BStask), max(Y_BStask), 100);
    xO = linspace(min(O_BStask), max(O_BStask), 100);

    plot(xY, polyval(pY_lin, xY), '-', 'Color', cols(1,:), 'LineWidth', 1.5)
    plot(xO, polyval(pO_lin, xO), '-', 'Color', cols(2,:), 'LineWidth', 1.5)

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

    % Fit linear and quadratic models on raw data
    mdl_lin  = fitlm(BStask, y_all_i);
    mdl_quad = fitlm(BStask, y_all_i, 'quadratic');

    % Manual nested F-test: does quadratic improve over linear?
    SSE_lin  = mdl_lin.SSE;
    SSE_quad = mdl_quad.SSE;

    df_lin   = mdl_lin.DFE;   % residual df
    df_quad  = mdl_quad.DFE;

    df1      = df_lin - df_quad;   % added parameter(s), usually 1
    df2      = df_quad;

    F_quad   = ((SSE_lin - SSE_quad) / df1) / (SSE_quad / df2);
    F_quad   = max(F_quad, 0); % guard against tiny negative values
    p_quad   = 1 - fcdf(F_quad, df1, df2);

    R2_lin   = mdl_lin.Rsquared.Adjusted;
    R2_quad  = mdl_quad.Rsquared.Adjusted;
    dR2      = R2_quad - R2_lin;

    if useBinnedAll
        % Bin x and average y within bins
        edges = linspace(min(BStask), max(BStask), nBins+1);

        yBinMean = nan(1,nBins);
        yBinSEM  = nan(1,nBins);
        xBinMean = nan(1,nBins);
        nPerBin  = nan(1,nBins);

        binIdx = discretize(BStask, edges);

        % Make sure max value goes into last bin
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

        % Plot binned means
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

        % Optional smooth line through binned means for visualization
        if sum(valid) >= 4
            if i == 2
                pBin = polyfit(xBinMean(valid), yBinMean(valid), 2);
            else
                pBin = polyfit(xBinMean(valid), yBinMean(valid), 1);
            end
            yFitBin = polyval(pBin, xFit);
            plot(xFit, yFitBin, 'k--', 'LineWidth', 1.4)
        end

    else
        % Plot model chosen for each panel
        if i == 2
            yFit = predict(mdl_quad, xFit');
        else
            yFit = predict(mdl_lin, xFit');
        end
        plot(xFit, yFit, 'k', 'LineWidth', 1.8)
    end

    % Print stats to command window
    fprintf('\n=== %s ===\n', regime_labels{i});
    fprintf('r_YA = %.3f, p = %.4f\n', rY, pY);
    fprintf('r_OA = %.3f, p = %.4f\n', rO, pO);
    fprintf('r_all = %.3f, p = %.4f\n', r_all, p_all);
    fprintf('Linear adj R2 = %.3f\n', R2_lin);
    fprintf('Quadratic adj R2 = %.3f\n', R2_quad);
    fprintf('Delta R2 = %.3f\n', dR2);
    fprintf('Quadratic vs linear: F(%d,%d) = %.3f, p = %.4f\n', ...
        df1, df2, F_quad, p_quad);

    % -------------------------
    % Title + annotation
    % -------------------------
    title(regime_labels{i})

    if p_all < 0.001
        p_str = 'p < 0.001';
    else
        p_str = sprintf('p = %.3f', p_all);
    end

    corr_str = sprintf('r = %.2f\n%s', r_all, p_str);

    if p_quad < 0.001
        p_quad_str = 'p < 0.001';
    else
        p_quad_str = sprintf('p = %.3f', p_quad);
    end

    quad_str = sprintf('quad > lin:\nF(%d,%d)=%.2f\n%s, dR^2=%.03f', ...
        df1, df2, F_quad, p_quad_str, dR2);

    if i == 2
        txt = sprintf('%s\n\n%s', corr_str, quad_str);
    else
        txt = corr_str;
    end

    text(0.92, 0.08, txt, ...
      'Units','normalized', ...
      'HorizontalAlignment','right', ...
      'VerticalAlignment','bottom', ...
      'FontSize',9, ...
      'BackgroundColor','none', ...
      'Margin',2);

    xlabel('Brain score (z)')
    ylabel('Behavior')

    xlim(xlim_all)
    ylim(ylim_all)
    box off
    axis padded
end

% outfile = sprintf('LV%d_BehavPLS', LVsel);
% exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), ...
%     'ContentType', 'vector', 'BackgroundColor','white')