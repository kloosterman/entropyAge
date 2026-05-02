%% =========================================================
% Correlation gradient: stable -> flexible cognitive domains
% versus mean-centered PLS brain scores
%
% Requires:
%   stat      : struct with stat.brainscores
%   clean_tbl : table with columns
%       Crystallized, Attention, LearningMemory,
%       WorkingMemory, FluidIntelligence, AgeGroup
%
% Notes:
% - Domain order is Stable -> Flexible
% - Plot shows correlations with bootstrap 95% CIs
% - Formal linear trend is tested across ordered domains
% ==========================================================

%% Settings
domainNames  = {'Crystallized','Attention','LearningMemory','WorkingMemory','FluidIntelligence'};
domainLabels = {'Crystallized','Attention','Learning & Memory','Working Memory','Fluid Intelligence'};

% Custom x positions so the 3 regimes take roughly equal visual space
% Stable       : Crystallized
% Balanced     : Attention, LearningMemory
% Flexible     : WorkingMemory, FluidIntelligence
% xPlot = [1.0 2.5 3.5 5.5 6.5];
xPlot = [1.0 2.5 3.5 4.5 5.5];

nBoot = 5000;
alpha = 0.05;
rng(1); % reproducible bootstrap

% Colors
col_all = [0 0 0];
col_Y   = [1.00 0.60 0.60];
col_O   = [0.60 0.75 0.95];

%% Extract brain scores
Y_BS   = stat.brainscores{1}(:,1);
O_BS   = stat.brainscores{2}(:,1);
ALL_BS = [Y_BS; O_BS];

%% Extract group indices
isY = clean_tbl.AgeGroup == "Young";
isO = clean_tbl.AgeGroup == "Older";

% Sanity check
assert(sum(isY) == numel(Y_BS), 'Young group size in clean_tbl does not match stat.brainscores{1}.');
assert(sum(isO) == numel(O_BS), 'Older group size in clean_tbl does not match stat.brainscores{2}.');

%% Preallocate
nDom = numel(domainNames);

r_all = nan(1,nDom);
r_Y   = nan(1,nDom);
r_O   = nan(1,nDom);

p_all = nan(1,nDom);
p_Y   = nan(1,nDom);
p_O   = nan(1,nDom);

ci_all = nan(2,nDom);
ci_Y   = nan(2,nDom);
ci_O   = nan(2,nDom);

%% Compute correlations and bootstrap CIs
for i = 1:nDom

    X_all = clean_tbl.(domainNames{i});
    X_Y   = X_all(isY);
    X_O   = X_all(isO);

    % Observed correlations
    [r_all(i), p_all(i)] = corr(ALL_BS, X_all, 'rows','pairwise', 'type','Pearson');
    [r_Y(i),   p_Y(i)]   = corr(Y_BS,   X_Y,   'rows','pairwise', 'type','Pearson');
    [r_O(i),   p_O(i)]   = corr(O_BS,   X_O,   'rows','pairwise', 'type','Pearson');

    % ----- Bootstrap ALL
    valid = ~isnan(ALL_BS) & ~isnan(X_all);
    x = X_all(valid);
    y = ALL_BS(valid);
    n = numel(x);

    boot_r = nan(nBoot,1);
    for b = 1:nBoot
        idx = randi(n,n,1);
        boot_r(b) = corr(y(idx), x(idx), 'type','Pearson');
    end
    ci_all(:,i) = prctile(boot_r, [100*alpha/2, 100*(1-alpha/2)]);

    % ----- Bootstrap YOUNG
    valid = ~isnan(Y_BS) & ~isnan(X_Y);
    x = X_Y(valid);
    y = Y_BS(valid);
    n = numel(x);

    boot_r = nan(nBoot,1);
    for b = 1:nBoot
        idx = randi(n,n,1);
        boot_r(b) = corr(y(idx), x(idx), 'type','Pearson');
    end
    ci_Y(:,i) = prctile(boot_r, [100*alpha/2, 100*(1-alpha/2)]);

    % ----- Bootstrap OLDER
    valid = ~isnan(O_BS) & ~isnan(X_O);
    x = X_O(valid);
    y = O_BS(valid);
    n = numel(x);

    boot_r = nan(nBoot,1);
    for b = 1:nBoot
        idx = randi(n,n,1);
        boot_r(b) = corr(y(idx), x(idx), 'type','Pearson');
    end
    ci_O(:,i) = prctile(boot_r, [100*alpha/2, 100*(1-alpha/2)]);
end

%% Convert CIs to asymmetric error bars
err_all_low = r_all - ci_all(1,:);
err_all_up  = ci_all(2,:) - r_all;

err_Y_low = r_Y - ci_Y(1,:);
err_Y_up  = ci_Y(2,:) - r_Y;

err_O_low = r_O - ci_O(1,:);
err_O_up  = ci_O(2,:) - r_O;

%% =========================================================
% Formal linear trend test across ordered domains
%
% Tests whether correlation strength changes linearly from
% stable -> flexible using the 5 ordered domains.
%
% NOTE:
% With only 5 domains, this is a descriptive/regime-level test.
% It treats the 5 domain-wise correlations as the observations.
% ==========================================================

ord = (1:nDom)';

% All subjects
mdl_all = fitlm(ord, r_all');
b_all   = mdl_all.Coefficients.Estimate(2);
t_all   = mdl_all.Coefficients.tStat(2);
pLin_all = mdl_all.Coefficients.pValue(2);
R2_all  = mdl_all.Rsquared.Ordinary;

% Young
mdl_Y = fitlm(ord, r_Y');
b_Y   = mdl_Y.Coefficients.Estimate(2);
t_Y   = mdl_Y.Coefficients.tStat(2);
pLin_Y = mdl_Y.Coefficients.pValue(2);
R2_Y  = mdl_Y.Rsquared.Ordinary;

% Older
mdl_O = fitlm(ord, r_O');
b_O   = mdl_O.Coefficients.Estimate(2);
t_O   = mdl_O.Coefficients.tStat(2);
pLin_O = mdl_O.Coefficients.pValue(2);
R2_O  = mdl_O.Rsquared.Ordinary;

%% Optional: test whether slopes differ between Young and Older
% Stack 10 observations: 5 domains x 2 groups
r_stack     = [r_Y(:); r_O(:)];
ord_stack   = repmat(ord, 2, 1);
group_stack = categorical([repmat("Young",nDom,1); repmat("Older",nDom,1)]);

tbl_trend = table(r_stack, ord_stack, group_stack, ...
    'VariableNames', {'r','Order','Group'});

mdl_inter = fitlm(tbl_trend, 'r ~ Order * Group');

% Interaction term = difference in slope between groups
coefNames = mdl_inter.CoefficientNames;
idx_inter = find(strcmp(coefNames, 'Order:Group_Older'));

if ~isempty(idx_inter)
    b_inter    = mdl_inter.Coefficients.Estimate(idx_inter);
    t_inter    = mdl_inter.Coefficients.tStat(idx_inter);
    p_inter    = mdl_inter.Coefficients.pValue(idx_inter);
else
    b_inter = NaN; t_inter = NaN; p_inter = NaN;
end

%% Print stats
fprintf('\n=== Domain-wise correlations with mcPLS brain score ===\n');
for i = 1:nDom
    fprintf('%s:\n', domainLabels{i});
    fprintf('  ALL   r = %+0.3f, p = %0.4f, 95%% CI [%+0.3f, %+0.3f]\n', ...
        r_all(i), p_all(i), ci_all(1,i), ci_all(2,i));
    fprintf('  YOUNG r = %+0.3f, p = %0.4f, 95%% CI [%+0.3f, %+0.3f]\n', ...
        r_Y(i), p_Y(i), ci_Y(1,i), ci_Y(2,i));
    fprintf('  OLDER r = %+0.3f, p = %0.4f, 95%% CI [%+0.3f, %+0.3f]\n', ...
        r_O(i), p_O(i), ci_O(1,i), ci_O(2,i));
end

fprintf('\n=== Linear trend across domains (Stable -> Flexible) ===\n');
fprintf('ALL   : slope = %+0.3f, t(3) = %+0.3f, p = %0.4f, R^2 = %0.3f\n', ...
    b_all, t_all, pLin_all, R2_all);
fprintf('YOUNG : slope = %+0.3f, t(3) = %+0.3f, p = %0.4f, R^2 = %0.3f\n', ...
    b_Y, t_Y, pLin_Y, R2_Y);
fprintf('OLDER : slope = %+0.3f, t(3) = %+0.3f, p = %0.4f, R^2 = %0.3f\n', ...
    b_O, t_O, pLin_O, R2_O);

fprintf('\n=== Group difference in linear trend ===\n');
fprintf('Order x Group interaction: b = %+0.3f, t = %+0.3f, p = %0.4f\n', ...
    b_inter, t_inter, p_inter);

%% =========================================================
% Plot
% =========================================================
% figure('Units','centimeters','Position',[5 5 12 7]); hold on
nexttile(7,[1 6]); 
hold on

% Young
errorbar(xPlot, r_Y, err_Y_low, err_Y_up, 'o-', ...
    'Color', col_Y, ...
    'MarkerFaceColor', col_Y, ...
    'MarkerEdgeColor', col_Y, ...
    'LineWidth', 0.5, ...
    'MarkerSize', 4, ...
    'CapSize', 0);

% Older
errorbar(xPlot, r_O, err_O_low, err_O_up, 'o-', ...
    'Color', col_O, ...
    'MarkerFaceColor', col_O, ...
    'MarkerEdgeColor', col_O, ...
    'LineWidth', 0.5, ...
    'MarkerSize', 4, ...
    'CapSize', 0);

% All subjects
errorbar(xPlot, r_all, err_all_low, err_all_up, 'o-', ...
    'Color', col_all, ...
    'MarkerFaceColor', col_all, ...
    'MarkerEdgeColor', col_all, ...
    'LineWidth', 1.25, ...
    'MarkerSize', 4, ...
    'CapSize', 0);

% Zero line
yline(0,'k:','LineWidth',0.5);

% Formatting
xlim([0.5 6.0]);
xticks(xPlot);
xticklabels(domainLabels);
xtickangle(25);

ylabel('Correlation with mcPLS brain score');

box off
ax = gca;
ax.FontSize = 8;
ax.TickDir = 'out';
ax.LineWidth = 0.75;

% Regime labels above the thirds
ylim([-1 1]);
yl = ylim;
yr = yl(2) - yl(1);

text(1.0, yl(2) + 0.01*yr, 'Stable', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',8);

text(3.0, yl(2) + 0.01*yr, 'Balanced', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',8);

text(5.0, yl(2) + 0.01*yr, 'Flexible', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontWeight','bold', ...
    'FontSize',8);

% ylim([yl(1), yl(2) + 0.16*yr]);

% Optional overall axis label above
% text(3.5, yl(2) + 0.15*yr, 'Cognitive processing gradient', ...
%     'HorizontalAlignment','center', ...
%     'VerticalAlignment','bottom', ...
%     'FontSize',8);

% Legend
lgd = legend({'Young','Older', 'All'}, ...
    'Location','northwest', ...
    'Box','off');
lgd.TextColor = 'k';

%% Optional title
% title('Stable-to-flexible gradient in brain-behavior coupling');