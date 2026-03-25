folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

% load stat_mse_2group_blink.mat
% load stat_mse_2group_blink_modulation.mat
load colormap_jetlightgray.mat

%% =========================
% SETTINGS
% =========================

% Choose which behavioral variables to include
% Indices refer to ORIGINAL behavioral variable order:
% 1 d2
% 2 VLMT 1–5
% 3 VLMT 1
% 4 VLMT 5
% 5 VLMT Interf.
% 6 VLMT Delayed
% 7 WMT-2
% 8 Digit span total
% 9 Digit span forw.
% 10 Digit span backw.
% 11 MWT-B
corrtype = "Pearson"; % pearson Spearman
plotinfig = 1

plotsel = 1
% Pretty labels in ORIGINAL order
pretty_labels_all = {
  'd2'
  'VLMT 1–5'
  'VLMT 1'
  'VLMT 5'
  'VLMT Interf.'
  'VLMT Delayed'
  'WMT-2'
  'Digit span'
  'Digit span forw.'
  'Digit span backw.'
  'MWT-B'
  };
if plotsel
  behavoi = [8 2 11];   % e.g. WM, Learning & Memory, Crystallized IQ
  domain_labels = {
    'Working Memory'
    'Learning & Memory'
    'Crystallized IQ'
    };
else
  behavoi = 1:11;   % e.g. WM, Learning & Memory, Crystallized IQ
  domain_labels = pretty_labels_all;
end

outfile = 'brain_behav_agediff_blink_oi';

%% =========================
% DATA
% =========================

behav_YA_all = stat_mse_2group_behav.results.stacked_behavdata(1:20,:);
behav_OA_all = stat_mse_2group_behav.results.stacked_behavdata(21:end,:);

brainscore_YA = stat_mse_2group_behav.brainscores{1}(:,1);
brainscore_OA = stat_mse_2group_behav.brainscores{2}(:,1);

% Select only variables of interest
behav_YA = behav_YA_all(:, behavoi);
behav_OA = behav_OA_all(:, behavoi);
pretty_labels = pretty_labels_all(behavoi);

nBehav = numel(behavoi);

%% =========================
% CORRELATIONS
% =========================

r_YA = nan(1,nBehav);
r_OA = nan(1,nBehav);

for ib = 1:nBehav
    r_YA(ib) = corr(behav_YA(:,ib), brainscore_YA, 'rows', 'complete', 'Type',corrtype);
    r_OA(ib) = corr(behav_OA(:,ib), brainscore_OA, 'rows', 'complete', 'Type',corrtype);
end

r_diff = r_OA - r_YA;

%% =========================
% BOOTSTRAP CI FOR YA / OA
% =========================

ul = stat_mse_2group_behav.boot_res.ulcorr;
ll = stat_mse_2group_behav.boot_res.llcorr;

% In full output: rows 1:11 = YA, 12:22 = OA
ul_YA_all = ul(1:11,1)';
ll_YA_all = ll(1:11,1)';
ul_OA_all = ul(12:22,1)';
ll_OA_all = ll(12:22,1)';

% Select only variables of interest
ul_YA = ul_YA_all(behavoi);
ll_YA = ll_YA_all(behavoi);
ul_OA = ul_OA_all(behavoi);
ll_OA = ll_OA_all(behavoi);

%% =========================
% BOOTSTRAP CI FOR OA - YA DIFFERENCE
% =========================

rng(1)
nBoot = 5000;

ll_diff = nan(1,nBehav);
ul_diff = nan(1,nBehav);

for ib = 1:nBehav

    xY = behav_YA(:,ib);
    yY = brainscore_YA;

    xO = behav_OA(:,ib);
    yO = brainscore_OA;

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

        rYb = corr(xY(idxY), yY(idxY));
        rOb = corr(xO(idxO), yO(idxO));

        boot_diff(b) = rOb - rYb;
    end

    ll_diff(ib) = prctile(boot_diff, 2.5);
    ul_diff(ib) = prctile(boot_diff, 97.5);
end

%% =========================
% ERROR BARS
% =========================

errY_low  = r_YA - ll_YA;
errY_high = ul_YA - r_YA;

errO_low  = r_OA - ll_OA;
errO_high = ul_OA - r_OA;

errD_low  = r_diff - ll_diff;
errD_high = ul_diff - r_diff;

%% =========================
% FISHER Z TEST FOR DIFFERENCE
% =========================

nY = size(behav_YA,1);
nO = size(behav_OA,1);

z_stat = nan(1,nBehav);
p_diff = nan(1,nBehav);

for i = 1:nBehav
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

col_YA_bar = [0.9 0.45 0.45];
col_OA_bar = [0.45 0.65 0.9];

col_YA_err = [0.75 0.15 0.15];
col_OA_err = [0.15 0.35 0.75];

col_diff = [0.15 0.15 0.15];

%% =========================
% PLOT
% =========================
if plotinfig
  nexttile(4,[1 3])
else
  f = figure;
  set(f, 'Units', 'centimeters')
  set(f, 'Position', [5 5 13 7.5])
end
hold on

R = [r_YA; r_OA]';
b = bar(R, 'grouped', 'BarWidth', 0.65);

b(1).FaceColor = col_YA_bar;
b(1).EdgeColor = 'none';

b(2).FaceColor = col_OA_bar;
b(2).EdgeColor = 'none';

xYA = b(1).XEndPoints;
xOA = b(2).XEndPoints;
xDiff = (xYA + xOA) / 2;

%% Error bars for YA/OA
errorbar(xYA, r_YA, errY_low, errY_high, ...
    'Color', col_YA_err, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 6)

errorbar(xOA, r_OA, errO_low, errO_high, ...
    'Color', col_OA_err, ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 6)

% %% Difference points
% errorbar(xDiff, r_diff, errD_low, errD_high, ...
%     'o', ...
%     'Color', col_diff, ...
%     'MarkerFaceColor', col_diff, ...
%     'MarkerEdgeColor', 'w', ...
%     'MarkerSize', 7, ...
%     'LineWidth', 1, ...
%     'CapSize', 6)

%% Optional text above bars
yl = [-1 1.1];
y_domain = yl(2) + 0.06;

for i = 1:nBehav
    text(i, y_domain, domain_labels{i}, ...
        'HorizontalAlignment', 'center', ...
        'FontSize', 9)
end

%% Significance markers
for i = 1:nBehav
    if sig_labels(i) ~= ""
        text(xDiff(i), r_OA(i) + 0.05, sig_labels(i), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 11, ...
            'FontWeight', 'bold')
    end
end

%% Axes
yline(0, 'k')

ylim(yl)
xlim([0.5 nBehav + 0.5])

xticks(1:nBehav)
xticklabels(pretty_labels)
% xtickangle(35)

% ylabel(sprintf('%s''s (r) and ∆ r', corrtype))
ylabel(sprintf('%s''s r', corrtype))

box off
set(gca, 'FontSize', 9)

% %% Export
% exportgraphics(gcf, fullfile(plotfolder, [outfile '.pdf']), 'ContentType', 'vector')
% exportgraphics(gcf, fullfile(plotfolder, [outfile '.svg']), 'ContentType', 'vector')