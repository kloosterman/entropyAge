%% Test LV1 x LV3 brain-score interaction predicting behavior
% Main question:
%   Does coordination between LV1 and LV3 predict behavior, and does that
%   LV1 x LV3 effect differ between young and older adults?
%
% Primary model per behavior outcome:
%   Behavior ~ LV1 * LV3 * Group

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

%% Data

lv1_YA = stat_mse_2group_behav.brainscores{1}(:,1);
lv3_YA = stat_mse_2group_behav.brainscores{1}(:,3);

lv1_OA = stat_mse_2group_behav.brainscores{2}(:,1);
lv3_OA = stat_mse_2group_behav.brainscores{2}(:,3);

nYA = numel(lv1_YA);
nOA = numel(lv1_OA);

if ~isfield(stat_mse_2group_behav, 'results') || ...
        ~isfield(stat_mse_2group_behav.results, 'stacked_behavdata')
    error('stat_mse_2group_behav.results.stacked_behavdata is missing.');
end

behav_all = stat_mse_2group_behav.results.stacked_behavdata;

if size(behav_all, 1) ~= nYA + nOA
    error('Behavior rows (%d) do not match brain-score rows (%d).', ...
        size(behav_all, 1), nYA + nOA);
end

% Original 5-domain order in stacked_behavdata:
% 1 WorkingMemory
% 2 FluidIntelligence
% 3 LearningMemory
% 4 Attention
% 5 Crystallized
domain_labels = {'WorkingMemory','FluidIntelligence','LearningMemory','Attention','Crystallized'};

behav_domain = behav_all;
behav_regime = nan(size(behav_all, 1), 3);
behav_regime(:,1) = mean(behav_all(:,[1 2]), 2, 'omitnan'); % Flexible
behav_regime(:,2) = mean(behav_all(:,[3 4]), 2, 'omitnan'); % Balanced
behav_regime(:,3) = behav_all(:,5);                         % Stable
regime_labels = {'Flexible','Balanced','Stable'};

Ymat = [behav_domain behav_regime];
outcome_labels = [domain_labels regime_labels];
outcome_type = [repmat("Domain", 1, numel(domain_labels)), ...
    repmat("Regime", 1, numel(regime_labels))];

LV1 = [lv1_YA; lv1_OA];
LV3 = [lv3_YA; lv3_OA];
LV1 = (LV1 - mean(LV1, 'omitnan')) ./ std(LV1, 'omitnan');
LV3 = (LV3 - mean(LV3, 'omitnan')) ./ std(LV3, 'omitnan');
Group = categorical([repmat("YA", nYA, 1); repmat("OA", nOA, 1)], ["YA","OA"]);

%% Fit models

nOut = size(Ymat, 2);

results = table( ...
    strings(nOut,1), strings(nOut,1), ...
    nan(nOut,1), nan(nOut,1), nan(nOut,1), ...
    nan(nOut,1), nan(nOut,1), nan(nOut,1), ...
    nan(nOut,1), nan(nOut,1), nan(nOut,1), ...
    nan(nOut,1), nan(nOut,1), ...
    nan(nOut,1), nan(nOut,1), ...
    'VariableNames', { ...
        'Outcome','Type','N','AdjR2','P_Model', ...
        'B_LV1xLV3_YA','SE_LV1xLV3_YA','P_LV1xLV3_YA', ...
        'B_LV1xLV3_OA','SE_LV1xLV3_OA','P_LV1xLV3_OA', ...
        'B_LV1xLV3xGroup','P_LV1xLV3xGroup', ...
        'F_LV1xLV3xGroup','P_Anova_LV1xLV3xGroup'});

results.B_LV1_YA = nan(nOut,1);
results.SE_LV1_YA = nan(nOut,1);
results.P_LV1_YA = nan(nOut,1);
results.B_LV1_OA = nan(nOut,1);
results.SE_LV1_OA = nan(nOut,1);
results.P_LV1_OA = nan(nOut,1);

results.B_LV3_YA = nan(nOut,1);
results.SE_LV3_YA = nan(nOut,1);
results.P_LV3_YA = nan(nOut,1);
results.B_LV3_OA = nan(nOut,1);
results.SE_LV3_OA = nan(nOut,1);
results.P_LV3_OA = nan(nOut,1);

models = cell(nOut, 1);
models_YA = cell(nOut, 1);
models_OA = cell(nOut, 1);

fprintf('\nLV1 x LV3 interaction predicting behavior\n');
fprintf('Full model: Behavior ~ LV1 * LV3 * Group\n');

for iout = 1:nOut
    Y = Ymat(:,iout);
    Y = (Y - mean(Y, 'omitnan')) ./ std(Y, 'omitnan');
    valid = ~isnan(Y) & ~isnan(LV1) & ~isnan(LV3);

    tbl = table(Y(valid), LV1(valid), LV3(valid), Group(valid), ...
        'VariableNames', {'Behavior','LV1','LV3','Group'});

    mdl = fitlm(tbl, 'Behavior ~ LV1 * LV3 * Group');
    models{iout} = mdl;

    coef_names = string(mdl.CoefficientNames);
    idx_ya_int = coef_names == "LV1:LV3";
    idx_3way = contains(coef_names, "LV1") & ...
        contains(coef_names, "LV3") & ...
        contains(coef_names, "Group");

    b_ya = nan;
    se_ya = nan;
    p_ya = nan;
    if any(idx_ya_int)
        b_ya = mdl.Coefficients.Estimate(find(idx_ya_int, 1));
        se_ya = mdl.Coefficients.SE(find(idx_ya_int, 1));
        p_ya = mdl.Coefficients.pValue(find(idx_ya_int, 1));
    end

    b_3way = nan;
    p_3way = nan;
    if any(idx_3way)
        b_3way = mdl.Coefficients.Estimate(find(idx_3way, 1));
        p_3way = mdl.Coefficients.pValue(find(idx_3way, 1));
    end

    idxYA = tbl.Group == "YA";
    idxOA = tbl.Group == "OA";

    mdl_YA = fitlm(tbl(idxYA,:), 'Behavior ~ LV1 * LV3');
    mdl_OA = fitlm(tbl(idxOA,:), 'Behavior ~ LV1 * LV3');
    models_YA{iout} = mdl_YA;
    models_OA{iout} = mdl_OA;

    coef_YA = string(mdl_YA.CoefficientNames);
    coef_OA = string(mdl_OA.CoefficientNames);

    idx_YA_group = coef_YA == "LV1:LV3";
    idx_OA_group = coef_OA == "LV1:LV3";
    idx_YA_LV1 = coef_YA == "LV1";
    idx_OA_LV1 = coef_OA == "LV1";
    idx_YA_LV3 = coef_YA == "LV3";
    idx_OA_LV3 = coef_OA == "LV3";

    if any(idx_YA_group)
        b_ya = mdl_YA.Coefficients.Estimate(find(idx_YA_group, 1));
        se_ya = mdl_YA.Coefficients.SE(find(idx_YA_group, 1));
        p_ya = mdl_YA.Coefficients.pValue(find(idx_YA_group, 1));
    end

    b_oa = nan;
    se_oa = nan;
    p_oa = nan;
    if any(idx_OA_group)
        b_oa = mdl_OA.Coefficients.Estimate(find(idx_OA_group, 1));
        se_oa = mdl_OA.Coefficients.SE(find(idx_OA_group, 1));
        p_oa = mdl_OA.Coefficients.pValue(find(idx_OA_group, 1));
    end

    b_lv1_ya = nan; se_lv1_ya = nan; p_lv1_ya = nan;
    b_lv1_oa = nan; se_lv1_oa = nan; p_lv1_oa = nan;
    b_lv3_ya = nan; se_lv3_ya = nan; p_lv3_ya = nan;
    b_lv3_oa = nan; se_lv3_oa = nan; p_lv3_oa = nan;

    if any(idx_YA_LV1)
        b_lv1_ya = mdl_YA.Coefficients.Estimate(find(idx_YA_LV1, 1));
        se_lv1_ya = mdl_YA.Coefficients.SE(find(idx_YA_LV1, 1));
        p_lv1_ya = mdl_YA.Coefficients.pValue(find(idx_YA_LV1, 1));
    end

    if any(idx_OA_LV1)
        b_lv1_oa = mdl_OA.Coefficients.Estimate(find(idx_OA_LV1, 1));
        se_lv1_oa = mdl_OA.Coefficients.SE(find(idx_OA_LV1, 1));
        p_lv1_oa = mdl_OA.Coefficients.pValue(find(idx_OA_LV1, 1));
    end

    if any(idx_YA_LV3)
        b_lv3_ya = mdl_YA.Coefficients.Estimate(find(idx_YA_LV3, 1));
        se_lv3_ya = mdl_YA.Coefficients.SE(find(idx_YA_LV3, 1));
        p_lv3_ya = mdl_YA.Coefficients.pValue(find(idx_YA_LV3, 1));
    end

    if any(idx_OA_LV3)
        b_lv3_oa = mdl_OA.Coefficients.Estimate(find(idx_OA_LV3, 1));
        se_lv3_oa = mdl_OA.Coefficients.SE(find(idx_OA_LV3, 1));
        p_lv3_oa = mdl_OA.Coefficients.pValue(find(idx_OA_LV3, 1));
    end

    anova_tbl = anova(mdl);
    anova_names = string(anova_tbl.Properties.RowNames);
    idx_anova_3way = contains(anova_names, "LV1") & ...
        contains(anova_names, "LV3") & ...
        contains(anova_names, "Group");

    F_3way = nan;
    p_anova_3way = nan;
    if any(idx_anova_3way)
        F_3way = anova_tbl.F(idx_anova_3way);
        p_anova_3way = anova_tbl.pValue(idx_anova_3way);
    end

    results.Outcome(iout) = outcome_labels{iout};
    results.Type(iout) = outcome_type(iout);
    results.N(iout) = height(tbl);
    results.AdjR2(iout) = mdl.Rsquared.Adjusted;
    results.P_Model(iout) = mdl.ModelFitVsNullModel.Pvalue;
    results.B_LV1xLV3_YA(iout) = b_ya;
    results.SE_LV1xLV3_YA(iout) = se_ya;
    results.P_LV1xLV3_YA(iout) = p_ya;
    results.B_LV1xLV3_OA(iout) = b_oa;
    results.SE_LV1xLV3_OA(iout) = se_oa;
    results.P_LV1xLV3_OA(iout) = p_oa;
    results.B_LV1xLV3xGroup(iout) = b_3way;
    results.P_LV1xLV3xGroup(iout) = p_3way;
    results.F_LV1xLV3xGroup(iout) = F_3way;
    results.P_Anova_LV1xLV3xGroup(iout) = p_anova_3way;

    results.B_LV1_YA(iout) = b_lv1_ya;
    results.SE_LV1_YA(iout) = se_lv1_ya;
    results.P_LV1_YA(iout) = p_lv1_ya;
    results.B_LV1_OA(iout) = b_lv1_oa;
    results.SE_LV1_OA(iout) = se_lv1_oa;
    results.P_LV1_OA(iout) = p_lv1_oa;

    results.B_LV3_YA(iout) = b_lv3_ya;
    results.SE_LV3_YA(iout) = se_lv3_ya;
    results.P_LV3_YA(iout) = p_lv3_ya;
    results.B_LV3_OA(iout) = b_lv3_oa;
    results.SE_LV3_OA(iout) = se_lv3_oa;
    results.P_LV3_OA(iout) = p_lv3_oa;

    fprintf('\n%s (%s)\n', outcome_labels{iout}, char(outcome_type(iout)));
    fprintf('  YA LV1xLV3: b = %+0.3f, p = %0.4f\n', b_ya, p_ya);
    fprintf('  OA LV1xLV3: b = %+0.3f, p = %0.4f\n', b_oa, p_oa);
    fprintf('  Group difference in LV1xLV3: b = %+0.3f, p = %0.4f\n', b_3way, p_3way);
end

%% FDR-correct primary tests

results.P_FDR_LV1xLV3_YA = mafdr(results.P_LV1xLV3_YA, 'BHFDR', true);
results.P_FDR_LV1xLV3_OA = mafdr(results.P_LV1xLV3_OA, 'BHFDR', true);
results.P_FDR_LV1xLV3xGroup = mafdr(results.P_LV1xLV3xGroup, 'BHFDR', true);

disp(results)

%% Contrast-based linear trend across regimes
% Regime coding requested here:
%   Stable   = -1
%   Balanced =  0
%   Flexible = +1
%
% With this coding, a positive LV1 x LV3 x RegimeLinear term means the
% LV1 x LV3 behavior coefficient increases toward Flexible. A negative term
% means it increases toward Stable.

regime_values = [-1 0 1];
regime_names = {'Stable','Balanced','Flexible'};
regime_cols = [3 2 1]; % columns in behav_regime are Flexible, Balanced, Stable

Y_long = [];
LV1_long = [];
LV3_long = [];
RegimeLinear_long = [];
Group_long = categorical(strings(0,1), ["YA","OA"]);

for ireg = 1:numel(regime_values)
    Yreg = behav_regime(:,regime_cols(ireg));
    Yreg = (Yreg - mean(Yreg, 'omitnan')) ./ std(Yreg, 'omitnan');

    Y_long = [Y_long; Yreg];
    LV1_long = [LV1_long; LV1];
    LV3_long = [LV3_long; LV3];
    RegimeLinear_long = [RegimeLinear_long; repmat(regime_values(ireg), size(LV1))];
    Group_long = [Group_long; Group];
end

valid_long = ~isnan(Y_long) & ~isnan(LV1_long) & ...
    ~isnan(LV3_long) & ~isnan(RegimeLinear_long);

trend_tbl = table( ...
    Y_long(valid_long), ...
    LV1_long(valid_long), ...
    LV3_long(valid_long), ...
    RegimeLinear_long(valid_long), ...
    Group_long(valid_long), ...
    'VariableNames', {'Behavior','LV1','LV3','RegimeLinear','Group'});

trend_model = fitlm(trend_tbl, 'Behavior ~ LV1 * LV3 * RegimeLinear * Group');

coef_names = string(trend_model.CoefficientNames);

idx_overall_gradient = coef_names == "LV1:LV3:RegimeLinear";
idx_group_gradient = contains(coef_names, "LV1") & ...
    contains(coef_names, "LV3") & ...
    contains(coef_names, "RegimeLinear") & ...
    contains(coef_names, "Group");

b_gradient_YA = nan;
se_gradient_YA = nan;
p_gradient_YA = nan;
if any(idx_overall_gradient)
    b_gradient_YA = trend_model.Coefficients.Estimate(find(idx_overall_gradient, 1));
    se_gradient_YA = trend_model.Coefficients.SE(find(idx_overall_gradient, 1));
    p_gradient_YA = trend_model.Coefficients.pValue(find(idx_overall_gradient, 1));
end

b_gradient_groupdiff = nan;
se_gradient_groupdiff = nan;
p_gradient_groupdiff = nan;
if any(idx_group_gradient)
    b_gradient_groupdiff = trend_model.Coefficients.Estimate(find(idx_group_gradient, 1));
    se_gradient_groupdiff = trend_model.Coefficients.SE(find(idx_group_gradient, 1));
    p_gradient_groupdiff = trend_model.Coefficients.pValue(find(idx_group_gradient, 1));
end

trend_tbl_YA = trend_tbl(trend_tbl.Group == "YA", :);
trend_tbl_OA = trend_tbl(trend_tbl.Group == "OA", :);

trend_model_YA = fitlm(trend_tbl_YA, 'Behavior ~ LV1 * LV3 * RegimeLinear');
trend_model_OA = fitlm(trend_tbl_OA, 'Behavior ~ LV1 * LV3 * RegimeLinear');

coef_YA = string(trend_model_YA.CoefficientNames);
coef_OA = string(trend_model_OA.CoefficientNames);
idx_YA_gradient = coef_YA == "LV1:LV3:RegimeLinear";
idx_OA_gradient = coef_OA == "LV1:LV3:RegimeLinear";

if any(idx_YA_gradient)
    b_gradient_YA = trend_model_YA.Coefficients.Estimate(find(idx_YA_gradient, 1));
    se_gradient_YA = trend_model_YA.Coefficients.SE(find(idx_YA_gradient, 1));
    p_gradient_YA = trend_model_YA.Coefficients.pValue(find(idx_YA_gradient, 1));
end

b_gradient_OA = nan;
se_gradient_OA = nan;
p_gradient_OA = nan;
if any(idx_OA_gradient)
    b_gradient_OA = trend_model_OA.Coefficients.Estimate(find(idx_OA_gradient, 1));
    se_gradient_OA = trend_model_OA.Coefficients.SE(find(idx_OA_gradient, 1));
    p_gradient_OA = trend_model_OA.Coefficients.pValue(find(idx_OA_gradient, 1));
end

anova_trend = anova(trend_model);
anova_names = string(anova_trend.Properties.RowNames);
idx_anova_overall = contains(anova_names, "LV1") & ...
    contains(anova_names, "LV3") & ...
    contains(anova_names, "RegimeLinear") & ...
    ~contains(anova_names, "Group");
idx_anova_groupdiff = contains(anova_names, "LV1") & ...
    contains(anova_names, "LV3") & ...
    contains(anova_names, "RegimeLinear") & ...
    contains(anova_names, "Group");

F_gradient_overall = nan;
p_anova_gradient_overall = nan;
if any(idx_anova_overall)
    F_gradient_overall = anova_trend.F(find(idx_anova_overall, 1));
    p_anova_gradient_overall = anova_trend.pValue(find(idx_anova_overall, 1));
end

F_gradient_groupdiff = nan;
p_anova_gradient_groupdiff = nan;
if any(idx_anova_groupdiff)
    F_gradient_groupdiff = anova_trend.F(find(idx_anova_groupdiff, 1));
    p_anova_gradient_groupdiff = anova_trend.pValue(find(idx_anova_groupdiff, 1));
end

trend_results = table( ...
    "Stable_minus1_Balanced0_Flexible_plus1", ...
    height(trend_tbl), ...
    trend_model.Rsquared.Adjusted, ...
    b_gradient_YA, se_gradient_YA, p_gradient_YA, ...
    b_gradient_OA, se_gradient_OA, p_gradient_OA, ...
    b_gradient_groupdiff, se_gradient_groupdiff, p_gradient_groupdiff, ...
    F_gradient_overall, p_anova_gradient_overall, ...
    F_gradient_groupdiff, p_anova_gradient_groupdiff, ...
    'VariableNames', { ...
        'RegimeCoding','N','AdjR2', ...
        'B_LV1xLV3xRegime_YA','SE_LV1xLV3xRegime_YA','P_LV1xLV3xRegime_YA', ...
        'B_LV1xLV3xRegime_OA','SE_LV1xLV3xRegime_OA','P_LV1xLV3xRegime_OA', ...
        'B_LV1xLV3xRegime_x_Group','SE_LV1xLV3xRegime_x_Group','P_LV1xLV3xRegime_x_Group', ...
        'F_LV1xLV3xRegime','P_Anova_LV1xLV3xRegime', ...
        'F_LV1xLV3xRegime_x_Group','P_Anova_LV1xLV3xRegime_x_Group'});

fprintf('\nContrast-based linear trend across regimes\n');
fprintf('Coding: Stable = -1, Balanced = 0, Flexible = +1\n');
fprintf('YA LV1xLV3 gradient: b = %+0.3f, p = %0.4f\n', ...
    b_gradient_YA, p_gradient_YA);
fprintf('OA LV1xLV3 gradient: b = %+0.3f, p = %0.4f\n', ...
    b_gradient_OA, p_gradient_OA);
fprintf('YA/OA gradient difference: b = %+0.3f, p = %0.4f\n', ...
    b_gradient_groupdiff, p_gradient_groupdiff);

disp(trend_results)

%% Summary plot: group-specific LV1 x LV3 coefficients

f = figure('Color', 'w');
f.Position = [200 200 760 360];

x = 1:nOut;
bar_width = 0.38;
x_YA_bar = x(:) - bar_width/2;
x_OA_bar = x(:) + bar_width/2;
hold on

hBarYA = bar(x - bar_width/2, results.B_LV1xLV3_YA, bar_width, ...
    'FaceColor', [0.85 0.25 0.25], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.75);
hBarOA = bar(x + bar_width/2, results.B_LV1xLV3_OA, bar_width, ...
    'FaceColor', [0.20 0.45 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.75);

errorbar(x_YA_bar, results.B_LV1xLV3_YA, results.SE_LV1xLV3_YA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);
errorbar(x_OA_bar, results.B_LV1xLV3_OA, results.SE_LV1xLV3_OA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);

yline(0, 'k-', 'LineWidth', 0.75);
xticks(x)
xticklabels(strrep(results.Outcome, '_', ' '))
xtickangle(35)
ylabel('LV1 x LV3 coefficient');
legend([hBarYA hBarOA], {'YA','OA'}, 'Location', 'best', 'Box', 'off');
box on
set(gca, 'FontSize', 9, 'TickDir', 'out', 'LineWidth', 1)
title('LV1 x LV3 interaction predicting behavior');

outfile = 'LV1xLV3_predict_behavior';

if exist('exportFigure', 'file')
    exportFigure(f, fullfile(plotfolder, [outfile '.pdf']), 'FontSize', 9);
    exportFigure(f, fullfile(plotfolder, [outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f, fullfile(plotfolder, [outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f, fullfile(plotfolder, [outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

%% Plot regime gradient coefficients

idx_regime = results.Type == "Regime";

f_trend = figure('Color', 'w');
f_trend.Position = [220 220 520 360];
hold on

x_reg = regime_values(:);
regime_results = results(idx_regime,:);
[~, plot_ord] = ismember(regime_names, cellstr(regime_results.Outcome));
y_YA = regime_results.B_LV1xLV3_YA(plot_ord);
y_OA = regime_results.B_LV1xLV3_OA(plot_ord);
se_YA = regime_results.SE_LV1xLV3_YA(plot_ord);
se_OA = regime_results.SE_LV1xLV3_OA(plot_ord);

hTrendYA = plot(x_reg, y_YA, 'o-', ...
    'Color', [0.85 0.25 0.25], ...
    'MarkerFaceColor', [0.85 0.25 0.25], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 7);
hTrendOA = plot(x_reg, y_OA, 'o-', ...
    'Color', [0.20 0.45 0.85], ...
    'MarkerFaceColor', [0.20 0.45 0.85], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 7);

errorbar(x_reg, y_YA, se_YA, ...
    'Color', [0.85 0.25 0.25], ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 5);
errorbar(x_reg, y_OA, se_OA, ...
    'Color', [0.20 0.45 0.85], ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 5);

yline(0, 'k-', 'LineWidth', 0.75);
xlim([-1.2 1.2])
xticks(regime_values)
xticklabels(regime_names)
ylabel('LV1 x LV3 coefficient');
legend([hTrendYA hTrendOA], {'YA','OA'}, 'Location', 'best', 'Box', 'off');
box on
set(gca, 'FontSize', 9, 'TickDir', 'out', 'LineWidth', 1)
title('Linear trend of LV1 x LV3 effect across regimes');

trend_outfile = 'LV1xLV3_regime_linear_trend';

if exist('exportFigure', 'file')
    exportFigure(f_trend, fullfile(plotfolder, [trend_outfile '.pdf']), 'FontSize', 9);
    exportFigure(f_trend, fullfile(plotfolder, [trend_outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f_trend, fullfile(plotfolder, [trend_outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f_trend, fullfile(plotfolder, [trend_outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

%% Separate plots for LV1 and LV3 main effects

f_lv1_bar = figure('Color', 'w');
f_lv1_bar.Position = [240 240 760 360];
hold on

hLV1BarYA = bar(x - bar_width/2, results.B_LV1_YA, bar_width, ...
    'FaceColor', [0.85 0.25 0.25], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.75);
hLV1BarOA = bar(x + bar_width/2, results.B_LV1_OA, bar_width, ...
    'FaceColor', [0.20 0.45 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.75);

errorbar(x_YA_bar, results.B_LV1_YA, results.SE_LV1_YA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);
errorbar(x_OA_bar, results.B_LV1_OA, results.SE_LV1_OA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);

yline(0, 'k-', 'LineWidth', 0.75);
xticks(x)
xticklabels(strrep(results.Outcome, '_', ' '))
xtickangle(35)
ylabel('LV1 coefficient');
legend([hLV1BarYA hLV1BarOA], {'YA','OA'}, 'Location', 'best', 'Box', 'off');
box on
set(gca, 'FontSize', 9, 'TickDir', 'out', 'LineWidth', 1)
title('LV1 main effect predicting behavior');

lv1_bar_outfile = 'LV1_main_effect_predict_behavior';

if exist('exportFigure', 'file')
    exportFigure(f_lv1_bar, fullfile(plotfolder, [lv1_bar_outfile '.pdf']), 'FontSize', 9);
    exportFigure(f_lv1_bar, fullfile(plotfolder, [lv1_bar_outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f_lv1_bar, fullfile(plotfolder, [lv1_bar_outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f_lv1_bar, fullfile(plotfolder, [lv1_bar_outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

f_lv3_bar = figure('Color', 'w');
f_lv3_bar.Position = [260 260 760 360];
hold on

hLV3BarYA = bar(x - bar_width/2, results.B_LV3_YA, bar_width, ...
    'FaceColor', [0.85 0.25 0.25], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.75);
hLV3BarOA = bar(x + bar_width/2, results.B_LV3_OA, bar_width, ...
    'FaceColor', [0.20 0.45 0.85], ...
    'EdgeColor', 'none', ...
    'FaceAlpha', 0.75);

errorbar(x_YA_bar, results.B_LV3_YA, results.SE_LV3_YA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);
errorbar(x_OA_bar, results.B_LV3_OA, results.SE_LV3_OA, ...
    'k', 'LineStyle', 'none', 'LineWidth', 0.8, 'CapSize', 4);

yline(0, 'k-', 'LineWidth', 0.75);
xticks(x)
xticklabels(strrep(results.Outcome, '_', ' '))
xtickangle(35)
ylabel('LV3 coefficient');
legend([hLV3BarYA hLV3BarOA], {'YA','OA'}, 'Location', 'best', 'Box', 'off');
box on
set(gca, 'FontSize', 9, 'TickDir', 'out', 'LineWidth', 1)
title('LV3 main effect predicting behavior');

lv3_bar_outfile = 'LV3_main_effect_predict_behavior';

if exist('exportFigure', 'file')
    exportFigure(f_lv3_bar, fullfile(plotfolder, [lv3_bar_outfile '.pdf']), 'FontSize', 9);
    exportFigure(f_lv3_bar, fullfile(plotfolder, [lv3_bar_outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f_lv3_bar, fullfile(plotfolder, [lv3_bar_outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f_lv3_bar, fullfile(plotfolder, [lv3_bar_outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

f_lv1_trend = figure('Color', 'w');
f_lv1_trend.Position = [280 280 520 360];
hold on

y_LV1_YA = regime_results.B_LV1_YA(plot_ord);
y_LV1_OA = regime_results.B_LV1_OA(plot_ord);
se_LV1_YA = regime_results.SE_LV1_YA(plot_ord);
se_LV1_OA = regime_results.SE_LV1_OA(plot_ord);

hLV1TrendYA = plot(x_reg, y_LV1_YA, 'o-', ...
    'Color', [0.85 0.25 0.25], ...
    'MarkerFaceColor', [0.85 0.25 0.25], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 7);
hLV1TrendOA = plot(x_reg, y_LV1_OA, 'o-', ...
    'Color', [0.20 0.45 0.85], ...
    'MarkerFaceColor', [0.20 0.45 0.85], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 7);

errorbar(x_reg, y_LV1_YA, se_LV1_YA, ...
    'Color', [0.85 0.25 0.25], ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 5);
errorbar(x_reg, y_LV1_OA, se_LV1_OA, ...
    'Color', [0.20 0.45 0.85], ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 5);

yline(0, 'k-', 'LineWidth', 0.75);
xlim([-1.2 1.2])
xticks(regime_values)
xticklabels(regime_names)
ylabel('LV1 coefficient');
legend([hLV1TrendYA hLV1TrendOA], {'YA','OA'}, 'Location', 'best', 'Box', 'off');
box on
set(gca, 'FontSize', 9, 'TickDir', 'out', 'LineWidth', 1)
title('LV1 main effect across regimes');

lv1_trend_outfile = 'LV1_main_effect_regime_trend';

if exist('exportFigure', 'file')
    exportFigure(f_lv1_trend, fullfile(plotfolder, [lv1_trend_outfile '.pdf']), 'FontSize', 9);
    exportFigure(f_lv1_trend, fullfile(plotfolder, [lv1_trend_outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f_lv1_trend, fullfile(plotfolder, [lv1_trend_outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f_lv1_trend, fullfile(plotfolder, [lv1_trend_outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

f_lv3_trend = figure('Color', 'w');
f_lv3_trend.Position = [300 300 520 360];
hold on

y_LV3_YA = regime_results.B_LV3_YA(plot_ord);
y_LV3_OA = regime_results.B_LV3_OA(plot_ord);
se_LV3_YA = regime_results.SE_LV3_YA(plot_ord);
se_LV3_OA = regime_results.SE_LV3_OA(plot_ord);

hLV3TrendYA = plot(x_reg, y_LV3_YA, 'o-', ...
    'Color', [0.85 0.25 0.25], ...
    'MarkerFaceColor', [0.85 0.25 0.25], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 7);
hLV3TrendOA = plot(x_reg, y_LV3_OA, 'o-', ...
    'Color', [0.20 0.45 0.85], ...
    'MarkerFaceColor', [0.20 0.45 0.85], ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 1.8, ...
    'MarkerSize', 7);

errorbar(x_reg, y_LV3_YA, se_LV3_YA, ...
    'Color', [0.85 0.25 0.25], ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 5);
errorbar(x_reg, y_LV3_OA, se_LV3_OA, ...
    'Color', [0.20 0.45 0.85], ...
    'LineStyle', 'none', ...
    'LineWidth', 1.0, ...
    'CapSize', 5);

yline(0, 'k-', 'LineWidth', 0.75);
xlim([-1.2 1.2])
xticks(regime_values)
xticklabels(regime_names)
ylabel('LV3 coefficient');
legend([hLV3TrendYA hLV3TrendOA], {'YA','OA'}, 'Location', 'best', 'Box', 'off');
box on
set(gca, 'FontSize', 9, 'TickDir', 'out', 'LineWidth', 1)
title('LV3 main effect across regimes');

lv3_trend_outfile = 'LV3_main_effect_regime_trend';

if exist('exportFigure', 'file')
    exportFigure(f_lv3_trend, fullfile(plotfolder, [lv3_trend_outfile '.pdf']), 'FontSize', 9);
    exportFigure(f_lv3_trend, fullfile(plotfolder, [lv3_trend_outfile '.png']), 'FontSize', 9, 'Resolution', 300);
else
    exportgraphics(f_lv3_trend, fullfile(plotfolder, [lv3_trend_outfile '.pdf']), ...
        'ContentType', 'vector', 'BackgroundColor', 'white');
    exportgraphics(f_lv3_trend, fullfile(plotfolder, [lv3_trend_outfile '.png']), ...
        'Resolution', 300, 'BackgroundColor', 'white');
end

writetable(results, fullfile(plotfolder, [outfile '_results.csv']));
writetable(trend_results, fullfile(plotfolder, [trend_outfile '_results.csv']));
save(fullfile(plotfolder, [outfile '_results.mat']), ...
    'results', 'trend_results', ...
    'models', 'models_YA', 'models_OA', ...
    'trend_model', 'trend_model_YA', 'trend_model_OA', ...
    'trend_tbl', 'anova_trend');

fprintf('\nSaved LV1 x LV3 behavior interaction results to %s\n', char(plotfolder));
