folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

%% Load behavioral data
behavfolder = '/Users/kloosterman/projectdata/EntropyAging/rsEEG-Daten/';
YA_tbl = readtable(fullfile(behavfolder, 'MA_Probanden_Neuropsychologie.xlsx'), 'Sheet', 'jüngere');
OA_tbl = readtable(fullfile(behavfolder, 'MA_Probanden_Neuropsychologie.xlsx'), 'Sheet', 'ältere');

load behav.mat

varNames = { ...
    'd2', ...
    'VLMT_1_5', ...
    'VLMT_1', ...
    'VLMT_5', ...
    'VLMT_Interf', ...
    'VLMT_Delayed', ...
    'WMT_2', ...
    'DigitSpan', ...
    'DigitSpan_Forw', ...
    'DigitSpan_Backw', ...
    'MWT_B'};

behav_tbl = array2table(behav, 'VariableNames', varNames);

nYA = 20;
nOA = 19;
n   = height(behav_tbl);

agegroup = strings(n,1);
agegroup(1:nYA) = "Young";
agegroup(nYA+1:end) = "Older";
behav_tbl.AgeGroup = categorical(agegroup);

%% Add subject indices (within group)

subjYA = (1:nYA)';          % 1..20
subjOA = (1:nOA)';          % 1..19

subjWithin = [subjYA; subjOA];

behav_tbl.SubjWithin = subjWithin;

% Optional: also keep a global index (useful for debugging)
behav_tbl.SubjGlobal = (1:height(behav_tbl))';

olind = behav_tbl.MWT_B <20;
behav_tbl.SubjWithin(olind)
behav_tbl.AgeGroup(olind)
%% -----------------------------------------
% Manual subject exclusion
% Indices are WITHIN group
% Example: oa_excl = [7] removes OA subject 7
%% -----------------------------------------
ya_excl = [];   % e.g. [3 8]
oa_excl = 5;   % 5 is the outlier in cryst

keepYA = true(nYA,1);
keepOA = true(nOA,1);

keepYA(ya_excl) = false;
keepOA(oa_excl) = false;

keepMask = [keepYA; keepOA];

behav_tbl = behav_tbl(keepMask,:);

fprintf('\n=== SUBJECT EXCLUSION SUMMARY ===\n');
fprintf('Young excluded (within-group indices): ');
if isempty(ya_excl)
    fprintf('none\n');
else
    fprintf('%d ', ya_excl);
    fprintf('\n');
end

fprintf('Older excluded (within-group indices): ');
if isempty(oa_excl)
    fprintf('none\n');
else
    fprintf('%d ', oa_excl);
    fprintf('\n');
end

fprintf('Remaining N: Young = %d, Older = %d, Total = %d\n', ...
    sum(keepYA), sum(keepOA), sum(keepMask));

%% Compute domain scores

% Working memory
WM_mat = [behav_tbl.DigitSpan_Forw, behav_tbl.DigitSpan_Backw];
WM_z = normalize(WM_mat);
behav_tbl.WorkingMemory = mean(WM_z, 2, 'omitnan');

% Fluid intelligence
behav_tbl.FluidIntelligence = behav_tbl.WMT_2;

% Learning & memory
VLMT_mat = [ ...
    behav_tbl.VLMT_1_5, ...
    behav_tbl.VLMT_1, ...
    behav_tbl.VLMT_5, ...
    behav_tbl.VLMT_Interf, ...
    behav_tbl.VLMT_Delayed];

VLMT_z = normalize(VLMT_mat);
behav_tbl.LearningMemory = mean(VLMT_z, 2, 'omitnan');

% Attention
behav_tbl.Attention = behav_tbl.d2;

% Crystallized intelligence
behav_tbl.Crystallized = behav_tbl.MWT_B;

%% Domain table in desired order (Stable → Flexible)
domain_tbl = behav_tbl(:, { ...
    'AgeGroup', ...
    'Crystallized', ...
    'Attention', ...
    'LearningMemory', ...
    'WorkingMemory', ...
    'FluidIntelligence'});

domainNames  = {'Crystallized','Attention','LearningMemory','WorkingMemory','FluidIntelligence'};

domainLabels = {'Crystallized','Attention','Learning & Memory','Working Memory','Fluid Intelligence'};

%% Z-score domains across all retained subjects
X = table2array(domain_tbl(:, 2:end));
Z = zscore(X);

clean_tbl = array2table(Z, 'VariableNames', domainNames);
clean_tbl.AgeGroup = domain_tbl.AgeGroup;

%% Group masks
isYoung = clean_tbl.AgeGroup == 'Young';
isOld   = clean_tbl.AgeGroup == 'Older';

%% Compute means, SEMs, and t-tests
meanY = zeros(1,5);
meanO = zeros(1,5);
semY  = zeros(1,5);
semO  = zeros(1,5);
pvals = zeros(1,5);

for i = 1:5
    dataY = clean_tbl{isYoung, domainNames{i}};
    dataO = clean_tbl{isOld,   domainNames{i}};

    meanY(i) = mean(dataY, 'omitnan');
    meanO(i) = mean(dataO, 'omitnan');

    semY(i) = std(dataY, 'omitnan') / sqrt(sum(~isnan(dataY)));
    semO(i) = std(dataO, 'omitnan') / sqrt(sum(~isnan(dataO)));

    [~, pvals(i)] = ttest2(dataY, dataO);
end

p_fdr = mafdr(pvals, 'BHFDR', true);

%% Brain score summary for appended BSV bars
BS_Y = stat_mse_2group_task.brainscores{1}(:,1);
BS_O = stat_mse_2group_task.brainscores{2}(:,1);

BS = zscore([BS_Y; BS_O]);
nBS_Y = numel(BS_Y);
BS_Y = BS(1:nBS_Y);
BS_O = BS(nBS_Y+1:end);

outBS_Y = abs(BS_Y) > 3;
outBS_O = abs(BS_O) > 3;
fprintf('Brain score outliers removed: Young=%d, Older=%d\n', sum(outBS_Y), sum(outBS_O));
BS_Y(outBS_Y) = NaN;
BS_O(outBS_O) = NaN;

meanBS_Y = mean(BS_Y, 'omitnan');
meanBS_O = mean(BS_O, 'omitnan');
semBS_Y = std(BS_Y, 'omitnan') / sqrt(sum(~isnan(BS_Y)));
semBS_O = std(BS_O, 'omitnan') / sqrt(sum(~isnan(BS_O)));
[~, pBS] = ttest2(BS_Y, BS_O);

%% Plot
nexttile(1,[1 3])
hold on

x = 1:6;
x_behav = 1:5;
x_bs = 6;
width = 0.35;
jitter = 0.08;

col_Y = [0.85 0.20 0.20];
col_O = [0.10 0.45 0.75];

meanY_plot = [meanY meanBS_Y];
meanO_plot = [meanO meanBS_O];
semY_plot = [semY semBS_Y];
semO_plot = [semO semBS_O];

% Bars
hY = bar(x - width/2, meanY_plot, width, ...
    'FaceColor', col_Y, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

hO = bar(x + width/2, meanO_plot, width, ...
    'FaceColor', col_O, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

% dividers between regimes
xline(1.5, '--', 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
xline(3.5, '--', 'LineWidth', 1, 'Color', [0.8 0.8 0.8])
xline(5.5, '-', 'LineWidth', 0.75, 'Color', [0.9 0.9 0.9])

% Individual points
for i = 1:5
    dataY_raw = clean_tbl{isYoung, domainNames{i}};
    dataO_raw = clean_tbl{isOld,   domainNames{i}};
    
    validY = ~isnan(dataY_raw);
    validO = ~isnan(dataO_raw);

    x_jit_Y = (x(i)-width/2) + (rand(sum(validY),1)-0.5)*2*jitter;
    x_jit_O = (x(i)+width/2) + (rand(sum(validO),1)-0.5)*2*jitter;

    scatter(x_jit_Y, dataY_raw(validY), ...
        4, col_Y, 'filled', ...
        'MarkerFaceAlpha', 0.4, ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.25);

    scatter(x_jit_O, dataO_raw(validO), ...
        4, col_O, 'filled', ...
        'MarkerFaceAlpha', 0.4, ...
        'MarkerEdgeColor', 'w', ...
        'LineWidth', 0.25);
end

x_jit_BS_Y = (x_bs-width/2) + (rand(sum(~isnan(BS_Y)),1)-0.5)*2*jitter;
x_jit_BS_O = (x_bs+width/2) + (rand(sum(~isnan(BS_O)),1)-0.5)*2*jitter;

scatter(x_jit_BS_Y, BS_Y(~isnan(BS_Y)), ...
    4, col_Y, 'filled', ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.25);

scatter(x_jit_BS_O, BS_O(~isnan(BS_O)), ...
    4, col_O, 'filled', ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.25);

% Error bars
errY = errorbar(x - width/2, meanY_plot, semY_plot, ...
    'Color', [0.35 0.35 0.35], ...
    'LineStyle', 'none', 'LineWidth', 0.5);
errO = errorbar(x + width/2, meanO_plot, semO_plot, ...
    'Color', [0.35 0.35 0.35], ...
    'LineStyle', 'none', 'LineWidth', 0.5);
set([errY errO], 'CapSize', 3);

% text(1.25, -2, 'Regimes:', 'FontSize', 7,'FontAngle','italic') %, 'Units','normalized'
% text(0.46, -2.7, 'Stable', 'FontSize', 7,'FontAngle','italic') %, 'Units','normalized'
% text(1.75, -2.7, 'Balanced', 'FontSize', 7,'FontAngle','italic') %, 'Units','normalized'
% text(3.77, -2.7, 'Flexible', 'FontSize', 7,'FontAngle','italic') %, 'Units','normalized'

%% Significance markers
sig_x = [1:5 x_bs];
sig_p = [p_fdr pBS];
sig_y = nan(size(sig_x));

for ii = 1:5
    dataY_raw = clean_tbl{isYoung, domainNames{ii}};
    dataO_raw = clean_tbl{isOld,   domainNames{ii}};
    sig_y(ii) = max([dataY_raw; dataO_raw], [], 'omitnan') + 0.22;
end
sig_y(end) = max([BS_Y; BS_O], [], 'omitnan') + 0.22;
sig_y(:) = sig_y(4); % shared height, matched to Working Memory

allVals = [clean_tbl{:, domainNames}(:); BS_Y(:); BS_O(:)];
yRange = max(allVals, [], 'omitnan') - min(allVals, [], 'omitnan');
if yRange == 0
    yRange = 1;
end

for ii = 1:numel(sig_x)
    if sig_p(ii) < 0.001
        stars = '***';
    elseif sig_p(ii) < 0.01
        stars = '**';
    elseif sig_p(ii) < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end

    plot([sig_x(ii)-width/2 sig_x(ii)+width/2], [sig_y(ii) sig_y(ii)], 'k', 'LineWidth', 1);
    text(sig_x(ii), sig_y(ii) + 0.03*yRange, stars, ...
        'HorizontalAlignment', 'center', 'FontSize', 9);
end

%% Formatting
xticks(x)
shortDomainLabels = {'Cryst.', 'Attention', 'Learning', 'Working', 'Fluid'};
xticklabels([shortDomainLabels {'Brain score'}])
xtickangle(25)

ylabel('Score (z)')
title('')

box off
ax = gca;
ax.XLim = [0.5 6.5];
if isempty(oa_excl)
  ax.YLim = [-3.4 max(sig_y) + 0.45];
else
  ax.YLim = [-3 max(sig_y) + 0.45];
end
ax.FontSize = 8;
ax.TickDir = 'out';
ax.YTick = [-2 0 2];
ax.XAxis.FontSize = 6;
panelA_YLim = ax.YLim;
panelA_YTick = ax.YTick;
panelA_Pos = ax.Position;

legend_y = -2.5;
plot(2.10, legend_y, 's', ...
    'MarkerFaceColor', col_Y, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', 4, ...
    'HandleVisibility', 'off');
text(2.22, legend_y, 'Young', ...
    'Color', 'k', ...
    'FontSize', 6, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle');

plot(4.25, legend_y, 's', ...
    'MarkerFaceColor', col_O, ...
    'MarkerEdgeColor', 'none', ...
    'MarkerSize', 4, ...
    'HandleVisibility', 'off');
text(4.37, legend_y, 'Older', ...
    'Color', 'k', ...
    'FontSize', 6, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle');
% exportgraphics(gcf, fullfile(plotfolder, 'domains_1col.pdf'), ...
%     'ContentType', 'vector', 'BackgroundColor', 'white');
% exportgraphics(gcf, fullfile(plotfolder, 'domains_1col.png'), ...
%     'BackgroundColor', 'white', 'Resolution', 600);
