%% plot correlations task PLS brainscores with raw behavior
% Figure
if not(plotinfig)
  f = figure;
  tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
end

% Data
Y_BStask = stat_mse_2group_task.brainscores{1}(:,1);
O_BStask = stat_mse_2group_task.brainscores{2}(:,1);
BStask   = [Y_BStask; O_BStask];

% behavoi = [7 2 11]; % Fluid intel, learning & memory, crystallized intelligence
% behaviorLabels = {'Fluid Intelligence', 'Learning & Memory', 'Crystallized Intelligence'};

behavoi = [8 2 11]; % Fluid intel, learning & memory, crystallized intelligence
 domain_labels = {
    'Working Memory'
    'Learning & Memory'
    'Crystallized IQ'
    };

% z-score each behavior variable across all subjects within test
behavdat = zscore(behav(:,behavoi));

nYA = length(Y_BStask);

% Identify OA rows in full dataset
idx_OA = (1:length(BStask)) > nYA;

% Work on full (z-scored) behavior
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

% Split behavior
behav_YA = behavdat(1:nYA,:);
behav_OA = behavdat(nYA+1:end,:);

% Common axis limits
x_all = BStask;
xlim_all = [min(x_all) max(x_all)];
ylim_all = [min(behavdat(:)) max(behavdat(:))];

cols = [1 0 0;          % Young = red
        0 0.447 0.741]; % Older = blue

for i = 1:3
    nexttile(i+6); hold on

    % Scatter
    h1 = scatter(Y_BStask, behav_YA(:,i), 30, 'filled', ...
        'MarkerFaceColor', cols(1,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.5);

    h2 = scatter(O_BStask, behav_OA(:,i), 30, 'filled', ...
        'MarkerFaceColor', cols(2,:), ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.5);

    % Linear fits within group
    pY = polyfit(Y_BStask, behav_YA(:,i), 1);
    pO = polyfit(O_BStask, behav_OA(:,i), 1);

    xY = linspace(min(Y_BStask), max(Y_BStask), 100);
    xO = linspace(min(O_BStask), max(O_BStask), 100);

    plot(xY, polyval(pY, xY), '-', 'Color', cols(1,:), 'LineWidth', 1.5)
    plot(xO, polyval(pO, xO), '-', 'Color', cols(2,:), 'LineWidth', 1.5)

    % Correlations
    [rY, pY] = corr(Y_BStask, behav_YA(:,i), 'rows', 'complete', 'Type',corrtype);
    [rO, pO] = corr(O_BStask, behav_OA(:,i), 'rows', 'complete', 'Type', corrtype);

    % Combined
    y_all_i = behavdat(:,i);
    [r_all, p_all] = corr(BStask, y_all_i, 'rows', 'complete');

    % Combined fit
    xFit = linspace(min(BStask), max(BStask), 200);

    if i == 2
        % quadratic for learning & memory
        pAll = polyfit(BStask, y_all_i, 2);
        yFit = polyval(pAll, xFit);
        plot(xFit, yFit, 'k', 'LineWidth', 1.8)
    else
        % linear for fluid + crystallized
        pAll = polyfit(BStask, y_all_i, 1);
        yFit = polyval(pAll, xFit);
        plot(xFit, yFit, 'k', 'LineWidth', 1.8)
    end

    % title(sprintf('%s\nr_{YA}=%.2f (p=%.3f), r_{OA}=%.2f (p=%.3f)\nr_{all}=%.2f (p=%.3f)', ...
    %   behaviorLabels{i}, rY, pY, rO, pO, r_all, p_all))
    % title(sprintf('%s\nr_{YA}=%.2f, r_{OA}=%.2f\nr_{all}=%.2f (p=%.3f)', ...
    %   domain_labels{i}, rY, rO, r_all, p_all))
    % title(sprintf('r_{YA}=%.2f, r_{OA}=%.2f\nr_{all}=%.2f (p=%.3f)', ...
    %   rY, rO, r_all, p_all))
    title(sprintf('r = %.2f, p=%.3f', ...
      r_all, p_all))

    xlabel('Brain score (z)')
    ylabel('Behavior (z)')

    xlim(xlim_all)
    ylim(ylim_all)
    box off
    axis padded

end

% legend([h1 h2], {'Young','Older'}, 'Location', 'best')

% test linear vs quadratic fit via f-test
% SSE and degrees of freedom
SSE_lin  = mdl_lin.SSE;
SSE_quad = mdl_quad.SSE;

df_lin  = mdl_lin.DFE;   % residual df
df_quad = mdl_quad.DFE;

% Difference in parameters
df1 = df_lin - df_quad;      % should be 1
df2 = df_quad;

% Partial F test
F = ((SSE_lin - SSE_quad)/df1) / (SSE_quad/df2);
p = 1 - fcdf(F, df1, df2);

fprintf('Linear vs quadratic: F(%d,%d) = %.3f, p = %.4f\n', df1, df2, F, p);
fprintf('Adj. R2 linear = %.3f\n', mdl_lin.Rsquared.Adjusted);
fprintf('Adj. R2 quadratic = %.3f\n', mdl_quad.Rsquared.Adjusted);

outfile = 'Main_figure';
exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), 'ContentType', 'vector')
% exportgraphics(gcf, fullfile(plotfolder, [outfile '.svg']), 'ContentType', 'vector')