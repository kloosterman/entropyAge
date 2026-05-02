folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

% load stat_mse_2group_blink.mat
% load stat_mse_2group_blink_modulation.mat
load colormap_jetlightgray.mat

%% =========================
% SETTINGS
% =========================

corrtype  = "Pearson";   % "Pearson" or "Spearman"
plotinfig = 1;
% LVsel     = 1;
doScatter = 1; % of BS vs behav

plot_mode = 'domain';    % 'regime' or 'domain'
outfile = ['brain_behav_agediff_blink_' plot_mode];

% Original 5-domain order in stacked_behavdata:
% 1 WorkingMemory
% 2 FluidIntelligence
% 3 LearningMemory
% 4 Attention
% 5 Crystallized

domain_labels = {'Working Memory','Fluid Intelligence','Learning & Memory','Attention','Crystallized'};
regime_labels = {'Stable','Balanced','Flexible'};

%% =========================
% DATA
% =========================

behav_YA_all = stat_mse_2group_behav.results.stacked_behavdata(1:20,:);
behav_OA_all = stat_mse_2group_behav.results.stacked_behavdata(21:end,:);

brainscore_YA = stat_mse_2group_behav.brainscores{1}(:,LVsel);
brainscore_OA = stat_mse_2group_behav.brainscores{2}(:,LVsel);

behav_ALL_all  = [behav_YA_all; behav_OA_all];
brainscore_ALL = [brainscore_YA; brainscore_OA];
group_ALL      = [zeros(size(brainscore_YA)); ones(size(brainscore_OA))];  % 0=YA, 1=OA

%% =========================
% SELECT MODE: DOMAIN OR REGIME
% =========================

switch lower(plot_mode)
  case 'domain'
    % Stable -> Balanced -> Flexible
    % Crystallized -> L&M -> Attention -> WM -> Fluid

    if LVsel == 1
      ord = [5 4 3 1 2];
    elseif LVsel == 3
      % also force LV3 into stable -> flexible ordering
      ord = [5 4 3 1 2];

    else
      ord = [5 3 4 1 2];
    end

    behav_YA  = behav_YA_all(:,ord);
    behav_OA  = behav_OA_all(:,ord);
    behav_ALL = behav_ALL_all(:,ord);
    xlabels   = domain_labels(ord);
    nBehav    = 5;

  case 'regime'
    % Stable, Balanced, Flexible

    behav_YA = nan(size(behav_YA_all,1), 3);
    behav_OA = nan(size(behav_OA_all,1), 3);
    behav_ALL = nan(size(behav_ALL_all,1), 3);

    behav_YA(:,1)  = behav_YA_all(:,5);                          % Stable
    behav_YA(:,2)  = mean(behav_YA_all(:,[3 4]), 2, 'omitnan');  % Balanced
    behav_YA(:,3)  = mean(behav_YA_all(:,[1 2]), 2, 'omitnan');  % Flexible

    behav_OA(:,1)  = behav_OA_all(:,5);
    behav_OA(:,2)  = mean(behav_OA_all(:,[3 4]), 2, 'omitnan');
    behav_OA(:,3)  = mean(behav_OA_all(:,[1 2]), 2, 'omitnan');

    behav_ALL(:,1) = behav_ALL_all(:,5);
    behav_ALL(:,2) = mean(behav_ALL_all(:,[3 4]), 2, 'omitnan');
    behav_ALL(:,3) = mean(behav_ALL_all(:,[1 2]), 2, 'omitnan');

    xlabels = regime_labels;
    nBehav  = 3;

  otherwise
    error('Unknown plot_mode. Use ''domain'' or ''regime''.')
end

% Ordered index for linear trend test (Stable -> Flexible)
ord_idx = (1:nBehav)';

%% =========================
% CORRELATIONS
% =========================

r_YA  = nan(1,nBehav);
r_OA  = nan(1,nBehav);
r_ALL = nan(1,nBehav);   % partial correlation controlling for group

for ib = 1:nBehav
  r_YA(ib) = corr(behav_YA(:,ib), brainscore_YA, 'rows','complete', 'Type', corrtype);
  r_OA(ib) = corr(behav_OA(:,ib), brainscore_OA, 'rows','complete', 'Type', corrtype);

  valid = ~isnan(behav_ALL(:,ib)) & ~isnan(brainscore_ALL) & ~isnan(group_ALL);
  r_ALL(ib) = partialcorr( ...
    behav_ALL(valid,ib), ...
    brainscore_ALL(valid), ...
    group_ALL(valid), ...
    'Type', corrtype);
end

r_diff = r_OA - r_YA;

%% =========================
% LINEAR TREND TEST
% =========================
% Tests whether the brain-behavior relationship becomes more positive
% from Stable -> Flexible, controlling for group

y_long      = [];
x_brain     = [];
x_domain    = [];
x_group     = [];

for ib = 1:nBehav
  y  = behav_ALL(:,ib);
  xB = brainscore_ALL;
  xD = repmat(ord_idx(ib), size(brainscore_ALL));
  xG = group_ALL;

  valid = ~isnan(y) & ~isnan(xB) & ~isnan(xD) & ~isnan(xG);

  y_long   = [y_long;   y(valid)];
  x_brain  = [x_brain;  xB(valid)];
  x_domain = [x_domain; xD(valid)];
  x_group  = [x_group;  xG(valid)];
end

% z-score continuous predictors for interpretable beta
x_brain_z  = zscore(x_brain);
x_domain_z = zscore(x_domain);

tblTrend = table(y_long, x_brain_z, x_domain_z, x_group, ...
  'VariableNames', {'Y','Brain','Domain','Group'});

mdl_trend = fitlm(tblTrend, 'Y ~ Brain*Domain + Group');

coefNames   = mdl_trend.CoefficientNames;
idxInteract = strcmp(coefNames, 'Brain:Domain');

beta_trend = mdl_trend.Coefficients.Estimate(idxInteract);
t_trend    = mdl_trend.Coefficients.tStat(idxInteract);
p_trend    = mdl_trend.Coefficients.pValue(idxInteract);

%% =========================
% 3-WAY INTERACTION TEST
% Brain × Domain × Group
% =========================

% Convert to categorical for proper ANOVA-style testing
DomainCat = categorical(x_domain);   % ordered domains
GroupCat  = categorical(x_group);    % 0=YA, 1=OA

% Build table
tbl3 = table(y_long, x_brain, DomainCat, GroupCat, ...
  'VariableNames', {'Y','Brain','Domain','Group'});

% Fit full interaction model
mdl_3way = fitlm(tbl3, 'Y ~ Brain*Domain*Group');

% ANOVA table (cleanest way to extract interaction test)
anova_tbl = anova(mdl_3way, 'summary');

% Display
disp(anova_tbl)

anova_tbl = anova(mdl_3way);

rowNames = anova_tbl.Properties.RowNames;

idx = strcmp(rowNames, 'Brain:Domain:Group');

F_3way = anova_tbl.F(idx);
p_3way = anova_tbl.pValue(idx);

fprintf('Brain × Domain × Group: F = %.2f, p = %.4f\n', F_3way, p_3way);

%% =========================
% CONTRAST: Stable vs Flexible
% =========================

% Identify domains
isStable   = x_domain == 1;           % after your reordering
isFlexible = x_domain == nBehav;

keep = isStable | isFlexible;

tbl_con = table( ...
  y_long(keep), ...
  x_brain(keep), ...
  x_group(keep), ...
  categorical(x_domain(keep)), ...
  'VariableNames', {'Y','Brain','Group','Domain'});

% Fit model
mdl_con = fitlm(tbl_con, 'Y ~ Brain*Domain*Group');

anova_con = anova(mdl_con);

disp(anova_con)
% disp(anova_con.Properties.RowNames)
%                  SumSq     DF    MeanSq       F         pValue
%                       ______    __    ______    _______    __________
%
% Brain                 24.798     1    24.798     2.6001       0.11135
% Group                 124.55     1    124.55      13.06    0.00056413
% Domain                7366.2     1    7366.2     772.38    1.5201e-39
% Brain:Group           48.807     1    48.807     5.1176      0.026791
% Brain:Domain          9.2377     1    9.2377    0.96861       0.32842
% Group:Domain           269.6     1     269.6     28.269    1.1978e-06
% Brain:Group:Domain    58.431     1    58.431     6.1267      0.015739
% Error                  667.6    70    9.5371
rowNames = anova_con.Properties.RowNames;

idx = contains(rowNames,'Brain') & contains(rowNames,'Domain') & contains(rowNames,'Group');

F_con = anova_con.F(idx);
p_con = anova_con.pValue(idx);

fprintf('Stable vs Flexible (Brain × Domain × Group): F = %.2f, p = %.4f\n', F_con, p_con);

%% =========================
% CONTRAST: Balanced vs Others
% =========================

isBalanced = x_domain == 2 | x_domain == 3;
isOther    = x_domain == 1 | x_domain == nBehav;

keep = isBalanced | isOther;

% recode domain: 0 = other, 1 = balanced
DomainBin = double(isBalanced(keep));

tbl_bal = table( ...
  y_long(keep), ...
  x_brain(keep), ...
  x_group(keep), ...
  DomainBin, ...
  'VariableNames', {'Y','Brain','Group','Balanced'});

tbl_bal.Balanced = categorical(tbl_bal.Balanced);

mdl_bal = fitlm(tbl_bal, 'Y ~ Brain*Balanced*Group');

anova_bal = anova(mdl_bal);

disp(anova_bal)

%% =========================
% CONFIDENCE INTERVALS
% =========================
% domain mode:
%   YA/OA CIs from PLS bootstrap output
%   ALL CI from direct bootstrap of partial correlation
%
% regime mode:
%   all CIs from direct bootstrap

rng(1)
nBoot = 5000;

ll_YA  = nan(1,nBehav);
ul_YA  = nan(1,nBehav);
ll_OA  = nan(1,nBehav);
ul_OA  = nan(1,nBehav);
ll_ALL = nan(1,nBehav);
ul_ALL = nan(1,nBehav);

if strcmpi(plot_mode, 'domain')
  % YA/OA from original PLS bootstrap output
  ul = stat_mse_2group_behav.boot_res.ulcorr;
  ll = stat_mse_2group_behav.boot_res.llcorr;

  ul_YA = ul(1:5,LVsel)';
  ll_YA = ll(1:5,LVsel)';
  ul_OA = ul(6:10,LVsel)';
  ll_OA = ll(6:10,LVsel)';

  ul_YA = ul_YA(ord);
  ll_YA = ll_YA(ord);
  ul_OA = ul_OA(ord);
  ll_OA = ll_OA(ord);

  % ALL-subject bootstrap CI for partial correlation controlling for group
  for ib = 1:nBehav
    xA = behav_ALL(:,ib);
    yA = brainscore_ALL;
    gA = group_ALL;

    validA = ~isnan(xA) & ~isnan(yA) & ~isnan(gA);
    xA = xA(validA);
    yA = yA(validA);
    gA = gA(validA);

    nA = numel(xA);
    boot_r = nan(nBoot,1);

    for b = 1:nBoot
      idxA = randi(nA, nA, 1);
      boot_r(b) = partialcorr(xA(idxA), yA(idxA), gA(idxA), 'Type', corrtype);
    end

    ll_ALL(ib) = prctile(boot_r, 2.5);
    ul_ALL(ib) = prctile(boot_r, 97.5);
  end

elseif strcmpi(plot_mode, 'regime')
  % direct bootstrap for YA / OA / ALL
  for ib = 1:nBehav

    % YA
    xY = behav_YA(:,ib);
    yY = brainscore_YA;
    validY = ~isnan(xY) & ~isnan(yY);
    xY = xY(validY);
    yY = yY(validY);
    nY = numel(xY);
    boot_rY = nan(nBoot,1);

    for b = 1:nBoot
      idx = randi(nY, nY, 1);
      boot_rY(b) = corr(xY(idx), yY(idx), 'Type', corrtype);
    end

    ll_YA(ib) = prctile(boot_rY, 2.5);
    ul_YA(ib) = prctile(boot_rY, 97.5);

    % OA
    xO = behav_OA(:,ib);
    yO = brainscore_OA;
    validO = ~isnan(xO) & ~isnan(yO);
    xO = xO(validO);
    yO = yO(validO);
    nO = numel(xO);
    boot_rO = nan(nBoot,1);

    for b = 1:nBoot
      idx = randi(nO, nO, 1);
      boot_rO(b) = corr(xO(idx), yO(idx), 'Type', corrtype);
    end

    ll_OA(ib) = prctile(boot_rO, 2.5);
    ul_OA(ib) = prctile(boot_rO, 97.5);

    % ALL partial correlation controlling for group
    xA = behav_ALL(:,ib);
    yA = brainscore_ALL;
    gA = group_ALL;

    validA = ~isnan(xA) & ~isnan(yA) & ~isnan(gA);
    xA = xA(validA);
    yA = yA(validA);
    gA = gA(validA);

    nA = numel(xA);
    boot_rA = nan(nBoot,1);

    for b = 1:nBoot
      idx = randi(nA, nA, 1);
      boot_rA(b) = partialcorr(xA(idx), yA(idx), gA(idx), 'Type', corrtype);
    end

    ll_ALL(ib) = prctile(boot_rA, 2.5);
    ul_ALL(ib) = prctile(boot_rA, 97.5);
  end
end

%% =========================
% BOOTSTRAP CI FOR OA - YA DIFFERENCE
% =========================

ll_diff = nan(1,nBehav);
ul_diff = nan(1,nBehav);

for ib = 1:nBehav
  xY = behav_YA(:,ib);  yY = brainscore_YA;
  xO = behav_OA(:,ib);  yO = brainscore_OA;

  validY = ~isnan(xY) & ~isnan(yY);
  validO = ~isnan(xO) & ~isnan(yO);

  xY = xY(validY); yY = yY(validY);
  xO = xO(validO); yO = yO(validO);

  nY = numel(xY);
  nO = numel(xO);

  boot_diff = nan(nBoot,1);

  for b = 1:nBoot
    idxY = randi(nY, nY, 1);
    idxO = randi(nO, nO, 1);

    rYb = corr(xY(idxY), yY(idxY), 'Type', corrtype);
    rOb = corr(xO(idxO), yO(idxO), 'Type', corrtype);

    boot_diff(b) = rOb - rYb;
  end

  ll_diff(ib) = prctile(boot_diff, 2.5);
  ul_diff(ib) = prctile(boot_diff, 97.5);
end

%% =========================
% ERROR BARS
% =========================

errY_low   = r_YA  - ll_YA;
errY_high  = ul_YA - r_YA;

errO_low   = r_OA  - ll_OA;
errO_high  = ul_OA - r_OA;

errA_low   = r_ALL  - ll_ALL;
errA_high  = ul_ALL - r_ALL;

%% =========================
% FISHER Z TEST FOR YA vs OA DIFFERENCE
% =========================

z_stat = nan(1,nBehav);
p_diff = nan(1,nBehav);

for i = 1:nBehav
  validY = ~isnan(behav_YA(:,i)) & ~isnan(brainscore_YA);
  validO = ~isnan(behav_OA(:,i)) & ~isnan(brainscore_OA);

  nY = sum(validY);
  nO = sum(validO);

  zY = atanh(r_YA(i));
  zO = atanh(r_OA(i));

  se = sqrt(1/(nY-3) + 1/(nO-3));
  z_stat(i) = (zO - zY) / se;
  p_diff(i) = 2 * (1 - normcdf(abs(z_stat(i))));
end

p_fdr = mafdr(p_diff, 'BHFDR', true);

%% =========================
% SIGNIFICANCE LABELS
% =========================

sig_labels = strings(1,nBehav);

for i = 1:nBehav
  if p_fdr(i) < .001
    sig_labels(i) = "***";
  elseif p_fdr(i) < .01
    sig_labels(i) = "**";
  elseif p_fdr(i) < .05
    sig_labels(i) = "*";
  end
end

%% =========================
% COLORS
% =========================

col_YA_bar  = [0.9 0.45 0.45];
col_OA_bar  = [0.45 0.65 0.9];

col_YA_err  = [0.75 0.15 0.15];
col_OA_err  = [0.15 0.35 0.75];

%% =========================
% PLOT
% =========================

if plotinfig
  axBar = nexttile(3,[1 2]);
  cla(axBar)
  hold(axBar,'on')
else
  f = figure;
  set(f, 'Units', 'centimeters')
  set(f, 'Position', [5 5 9 7.5])
end
hold on

R = [r_YA; r_OA]';
b = bar(R, 'grouped', 'BarWidth', 0.8);

b(1).FaceColor = col_YA_bar;   b(1).EdgeColor = 'none';
b(2).FaceColor = col_OA_bar;   b(2).EdgeColor = 'none';

xYA  = b(1).XEndPoints;
xOA  = b(2).XEndPoints;
xMid = 1:nBehav;

errorbar(xYA, r_YA, errY_low, errY_high, ...
  'Color', col_YA_err, 'LineStyle', 'none', 'LineWidth', 0.5, 'CapSize', 6);

errorbar(xOA, r_OA, errO_low, errO_high, ...
  'Color', col_OA_err, 'LineStyle', 'none', 'LineWidth', 0.5, 'CapSize', 6);

hAll = errorbar(xMid, r_ALL, errA_low, errA_high, ...
  '-o', ...
  'Color', [0 0 0], ...
  'MarkerFaceColor', [0 0 0], ...
  'MarkerEdgeColor', [0 0 0], ...
  'LineWidth', 1.2, ...
  'MarkerSize', 4, ...
  'CapSize', 6);

%% Significance markers for YA vs OA difference
for i = 1:nBehav
  if sig_labels(i) ~= ""
    ysig = max([r_YA(i)+errY_high(i), r_OA(i)+errO_high(i), r_ALL(i)+errA_high(i)]) + 0.05;
    text(mean([xYA(i) xOA(i)]), ysig, sig_labels(i), ...
      'HorizontalAlignment', 'center', ...
      'FontSize', 11, ...
      'FontWeight', 'bold')
  end
end

% %% Trend text
% txt = sprintf('Linear trend: \\beta = %.2f, t = %.2f, p = %.3f', beta_trend, t_trend, p_trend);
% text(0.02, 0.98, txt, 'Units','normalized', ...
%     'HorizontalAlignment','left', ...
%     'VerticalAlignment','top', ...
%     'FontSize', 7);

%% Axes
yline(0, 'k')
ylim([-1 1])
xlim([0.5 nBehav + 0.5])

xticks(1:nBehav)
xticklabels(xlabels)
ylabel(sprintf('%s''s r', corrtype))
xlabel('Cognitive domain (Stable \rightarrow Flexible)')

lgd = legend([b(1) b(2) hAll], {'Young','Older','Partial r | group'}, 'Location','southeast',  'Orientation','horizontal');
lgd.ItemTokenSize = [5 10];
legend boxoff

box off
ax = gca;
ax.LineWidth = 0.5;
ax.TickDir = 'out';
ax.FontSize = 8;

outfile = sprintf('LV%d_BehavPLS_%s', LVsel, plot_mode);
disp(outfile)
exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), 'ContentType', 'vector', 'BackgroundColor','white')
exportgraphics(gcf, fullfile(plotfolder, [outfile '.png']), 'BackgroundColor','white')

%% doScatter
idx_fluid = find(contains(lower(behav_labels), 'workingmemory'));   % fluid (WMT-2)
idx_cryst = find(contains(lower(behav_labels), 'crystallized'));   % crystallized (MWT-B)
doQuadratic = 0;

if doScatter
  behav_labels = domain_tbl.Properties.VariableNames(2:end);
  figure;

  nCols = 5;
  nRows = ceil(nBehav / nCols);
  tiledlayout(nRows, nCols, 'TileSpacing','compact','Padding','compact');

  for ib = 1:nBehav
    nexttile; hold on;

    % -------------------------
    % Scatter YA
    % -------------------------
    scatter(brainscore_YA, behav_YA(:,ib), ...
      30, [1 0 0], 'filled', 'MarkerFaceAlpha', 0.6);

    % -------------------------
    % YA fit
    % -------------------------
    valid_YA = ~isnan(brainscore_YA) & ~isnan(behav_YA(:,ib));

    if sum(valid_YA) > 3
      x = brainscore_YA(valid_YA);
      y = behav_YA(valid_YA,ib);

      if doQuadratic && ismember(ib, idx_fluid)
        % quadratic fit for YA fluid
        p = polyfit(x, y, 2);
      else
        % linear fit otherwise
        p = polyfit(x, y, 1);
      end

      xfit = linspace(min(x), max(x), 100);
      yfit = polyval(p, xfit);

      plot(xfit, yfit, 'Color', [1 0 0], 'LineWidth', 1.5);
    end

    % -------------------------
    % Scatter OA
    % -------------------------
    scatter(brainscore_OA, behav_OA(:,ib), ...
      30, [0 0.447 0.741], 'filled', 'MarkerFaceAlpha', 0.6);

    % -------------------------
    % OA fit
    % -------------------------
    valid_OA = ~isnan(brainscore_OA) & ~isnan(behav_OA(:,ib));

    if sum(valid_OA) > 3
      x = brainscore_OA(valid_OA);
      y = behav_OA(valid_OA,ib);

      if doQuadratic && ismember(ib, idx_cryst)
        % quadratic fit for OA crystallized
        p = polyfit(x, y, 2);
      else
        % linear fit otherwise
        p = polyfit(x, y, 1);
      end

      xfit = linspace(min(x), max(x), 100);
      yfit = polyval(p, xfit);

      plot(xfit, yfit, 'Color', [0 0.447 0.741], 'LineWidth', 1.5);
    end

    % -------------------------
    % Labels
    % -------------------------
    title(sprintf('%s\nYA r=%.2f | OA r=%.2f', ...
      behav_labels{ib}, r_YA(ib), r_OA(ib)), ...
      'Interpreter','none', 'FontSize',9);

    box on;
  end
  legend({'YA','YA fit','OA','OA fit'}, 'Location','southoutside', ...
    'Orientation','horizontal');
end

outfile = sprintf('LV%d_BehavPLS_%s_scatter', LVsel, plot_mode);
disp(outfile)
exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), 'ContentType', 'vector', 'BackgroundColor','white')
exportgraphics(gcf, fullfile(plotfolder, [outfile '.png']), 'BackgroundColor','white')
