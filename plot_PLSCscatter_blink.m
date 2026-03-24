folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

% load stat_mse_2group_blink.mat
% load stat_mse_2group_blink_modulation.mat

load colormap_jetlightgray.mat
%% Brain–behavior correlations: Young vs Older

% YA rows = 1:20
% OA rows = 21:end
% CI rows: 1:11 = YA, 12:22 = OA

%% Data
behav_YA = stat_mse_2group.results.stacked_behavdata(1:20,:);
behav_OA = stat_mse_2group.results.stacked_behavdata(21:end,:);
% behav_YA = stat_mse_2group.results.stacked_behavdata(1:19,:);
% behav_OA = stat_mse_2group.results.stacked_behavdata(20:end,:);

brainscore_YA = stat_mse_2group.brainscores{1}(:,1);
brainscore_OA = stat_mse_2group.brainscores{2}(:,1);

%% Remove YA subject with lowest brain score everywhere
% brainscore_YA_orig = brainscore_YA;
% [brain_min_YA, idx_remove_YA] = min(brainscore_YA_orig);
% 
% fprintf('Removed YA subject %d (lowest brain score = %.2f)\n', ...
%     idx_remove_YA, brain_min_YA);
% 
% brainscore_YA(idx_remove_YA) = [];
% behav_YA(idx_remove_YA,:) = [];

nBehav = size(behav_YA,2);

%% Correlations
r_YA = nan(1,nBehav);
r_OA = nan(1,nBehav);

for ib = 1:nBehav
    r_YA(ib) = corr(behav_YA(:,ib), brainscore_YA,'rows','complete');
    r_OA(ib) = corr(behav_OA(:,ib), brainscore_OA,'rows','complete');
end

r_diff = r_OA - r_YA;

%% Bootstrap CI for YA/OA
ul = stat_mse_2group.boot_res.ulcorr;
ll = stat_mse_2group.boot_res.llcorr;

ul_YA = ul(1:11,1)';
ll_YA = ll(1:11,1)';

ul_OA = ul(12:22,1)';
ll_OA = ll(12:22,1)';

% remove same YA subject from bootstrap-based CIs not possible post hoc,
% so these are retained from original solution as descriptive bounds

%% Bootstrap CI for OA–YA difference
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
        idxY = randi(nY,nY,1);
        idxO = randi(nO,nO,1);

        rYb = corr(xY(idxY),yY(idxY));
        rOb = corr(xO(idxO),yO(idxO));

        boot_diff(b) = rOb - rYb;
    end

    ll_diff(ib) = prctile(boot_diff,2.5);
    ul_diff(ib) = prctile(boot_diff,97.5);

end

%% Convert CI bounds to errorbars
errY_low  = r_YA - ll_YA;
errY_high = ul_YA - r_YA;

errO_low  = r_OA - ll_OA;
errO_high = ul_OA - r_OA;

errD_low  = r_diff - ll_diff;
errD_high = ul_diff - r_diff;

%% Fisher z-test for correlation differences
nY = size(behav_YA,1);
nO = size(behav_OA,1);

z_stat = nan(1,nBehav);
p_diff = nan(1,nBehav);

for i = 1:nBehav

    zY = atanh(r_YA(i));
    zO = atanh(r_OA(i));

    se = sqrt(1/(nY-3) + 1/(nO-3));

    z_stat(i) = (zO - zY)/se;

    p_diff(i) = 2*(1-normcdf(abs(z_stat(i))));

end

%% FDR correction
p_fdr = mafdr(p_diff,'BHFDR',true);

%% Significance labels
sig_labels = strings(1,nBehav);

for i=1:nBehav

    if p_fdr(i) < .001
        sig_labels(i)="***";
    elseif p_fdr(i) < .01
        sig_labels(i)="**";
    elseif p_fdr(i) < .08
        sig_labels(i)="*";
    end

end

%% Labels
pretty_labels = {
'd2'
'VLMT 1–5'
'VLMT 1'
'VLMT 5'
'VLMT Interf.'
'VLMT Delayed'
'WMT-2'
'Digit span total'
'Digit span forw.'
'Digit span backw.'
'MWT-B'
};

%% Reorder: VLMT first, then d2, then the rest
ord = [2 3 6 4 5 1 10 8 9 7 11];

r_YA    = r_YA(ord);
r_OA    = r_OA(ord);
r_diff  = r_diff(ord);

ll_YA   = ll_YA(ord);
ul_YA   = ul_YA(ord);
ll_OA   = ll_OA(ord);
ul_OA   = ul_OA(ord);

ll_diff = ll_diff(ord);
ul_diff = ul_diff(ord);

errY_low  = errY_low(ord);
errY_high = errY_high(ord);
errO_low  = errO_low(ord);
errO_high = errO_high(ord);
errD_low  = errD_low(ord);
errD_high = errD_high(ord);

z_stat     = z_stat(ord);
p_diff     = p_diff(ord);
p_fdr      = p_fdr(ord);
sig_labels = sig_labels(ord);

pretty_labels = pretty_labels(ord);

%% Colors
col_YA = [0.85 0.25 0.25];
col_OA = [0.20 0.45 0.85];

%% ------------------------------------------------------------------------
%% GLOBAL scatter plots + model comparison table
%% Brain score (x) vs raw behavior (y)
%% ------------------------------------------------------------------------

% Representative tests (raw column indices)
sel_idx = [2 1 10 7 11];   % VLMT 1-5, d2, Digit span backw., WMT-2, MWT-B

domain_names = { ...
    'Verbal memory', ...
    'Attention', ...
    'Working memory', ...
    'Fluid IQ', ...
    'Crystallized IQ'};

test_labels = { ...
    'VLMT 1–5', ...
    'd2', ...
    'Digit span backw.', ...
    'WMT-2', ...
    'MWT-B'};

%% Figure
f = figure;
set(f,'Units','centimeters')
set(f,'Position',[5 5 13 9.5])

tiledlayout(2,3,'TileSpacing','compact','Padding','compact');

% Store handles once for legend
hY = gobjects(1,1);
hO = gobjects(1,1);
hLY = gobjects(1,1);
hLO = gobjects(1,1);
hQ  = gobjects(1,1);

%% Output table
results_table = table;

for i = 1:numel(sel_idx)

    ib = sel_idx(i);

    % ----------------------------
    % Data
    % ----------------------------
    xY = brainscore_YA;
    yY = behav_YA(:,ib);

    xO = brainscore_OA;
    yO = behav_OA(:,ib);

    validY = ~isnan(xY) & ~isnan(yY);
    validO = ~isnan(xO) & ~isnan(yO);

    xY = xY(validY); yY = yY(validY);
    xO = xO(validO); yO = yO(validO);

    x_all = [xY; xO];
    y_all = [yY; yO];
    age   = [zeros(numel(xY),1); ones(numel(xO),1)];   % YA = 0, OA = 1

    % ----------------------------
    % Correlations
    % ----------------------------
    rY = corr(xY,yY,'rows','complete');
    rO = corr(xO,yO,'rows','complete');

    % ----------------------------
    % Models
    % ----------------------------
    T = table(y_all, x_all, age, 'VariableNames', {'behav','brain','age'});
    T.brain2 = T.brain.^2;

    mdl_lin  = fitlm(T, 'behav ~ brain');
    mdl_int  = fitlm(T, 'behav ~ brain * age');
    mdl_quad = fitlm(T, 'behav ~ brain + brain2');
    mdl_full = fitlm(T, 'behav ~ brain * age + brain2');

    % p-values safely by coefficient name
    coefnames_int  = mdl_int.CoefficientNames;
    coefnames_quad = mdl_quad.CoefficientNames;

    p_interaction = NaN;
    p_quadratic   = NaN;

    idx_int = find(strcmp(coefnames_int,'brain:age'));
    if ~isempty(idx_int)
        p_interaction = mdl_int.Coefficients.pValue(idx_int);
    end

    idx_quad = find(strcmp(coefnames_quad,'brain2'));
    if ~isempty(idx_quad)
        p_quadratic = mdl_quad.Coefficients.pValue(idx_quad);
    end

    % AIC deltas: positive means first-named model is better
    dAIC_quad_vs_int  = mdl_int.ModelCriterion.AIC  - mdl_quad.ModelCriterion.AIC;
    dAIC_quad_vs_lin  = mdl_lin.ModelCriterion.AIC  - mdl_quad.ModelCriterion.AIC;
    dAIC_full_vs_quad = mdl_quad.ModelCriterion.AIC - mdl_full.ModelCriterion.AIC;

    % ----------------------------
    % Save in table
    % ----------------------------
    results_table.Test{i,1}          = test_labels{i};
    results_table.r_YA(i,1)          = rY;
    results_table.r_OA(i,1)          = rO;

    results_table.AIC_lin(i,1)       = mdl_lin.ModelCriterion.AIC;
    results_table.AIC_int(i,1)       = mdl_int.ModelCriterion.AIC;
    results_table.AIC_quad(i,1)      = mdl_quad.ModelCriterion.AIC;
    results_table.AIC_full(i,1)      = mdl_full.ModelCriterion.AIC;

    results_table.AdjR2_lin(i,1)     = mdl_lin.Rsquared.Adjusted;
    results_table.AdjR2_int(i,1)     = mdl_int.Rsquared.Adjusted;
    results_table.AdjR2_quad(i,1)    = mdl_quad.Rsquared.Adjusted;
    results_table.AdjR2_full(i,1)    = mdl_full.Rsquared.Adjusted;

    results_table.p_interaction(i,1) = p_interaction;
    results_table.p_quadratic(i,1)   = p_quadratic;

    results_table.dAIC_quad_vs_int(i,1)  = dAIC_quad_vs_int;
    results_table.dAIC_quad_vs_lin(i,1)  = dAIC_quad_vs_lin;
    results_table.dAIC_full_vs_quad(i,1) = dAIC_full_vs_quad;

    % ----------------------------
    % Plot
    % ----------------------------
    nexttile
    hold on

    % points
    hY = scatter(xY,yY,30,'filled', ...
        'MarkerFaceColor',col_YA, ...
        'MarkerEdgeColor','w', ...
        'LineWidth',0.5);

    hO = scatter(xO,yO,30,'filled', ...
        'MarkerFaceColor',col_OA, ...
        'MarkerEdgeColor','w', ...
        'LineWidth',0.5);

    % group-wise linear fits
    pY = polyfit(xY,yY,1);
    pO = polyfit(xO,yO,1);

    xxY = linspace(min(xY),max(xY),100);
    xxO = linspace(min(xO),max(xO),100);
    xxA = linspace(min(x_all),max(x_all),200);

    hLY = plot(xxY,polyval(pY,xxY),'-','Color',col_YA,'LineWidth',1.3);
    hLO = plot(xxO,polyval(pO,xxO),'-','Color',col_OA,'LineWidth',1.3);

    % global quadratic fit
    pQ = polyfit(x_all,y_all,2);
    hQ = plot(xxA,polyval(pQ,xxA),'k--','LineWidth',1.8);

    % zero line
    xline(0,'k:','LineWidth',0.75)

    % labels
    title({domain_names{i}, test_labels{i}}, ...
        'FontSize',9,'FontWeight','normal')

    xlabel('Brain score')
    ylabel(test_labels{i})
    xlim([500 1500])

    % annotation
    text(0.03,0.97, ...
        sprintf('r_Y = %.2f\nr_O = %.2f\n\\DeltaAIC Q-I = %.2f', ...
        rY, rO, dAIC_quad_vs_int), ...
        'Units','normalized', ...
        'HorizontalAlignment','left', ...
        'VerticalAlignment','top', ...
        'FontSize',8)

    box off
    set(gca,'FontSize',8)

end

%% Legend / info tile
axL = nexttile(6);
axis(axL,'off')

legend(axL,[hY hO hLY hLO hQ], ...
    {'Young','Older','YA linear','OA linear','Global quadratic'}, ...
    'Location','northwest','Box','off','FontSize',8)

text(0,0.58, ...
    '\DeltaAIC Q-I > 0 means quadratic fits better than interaction.', ...
    'Parent',axL,'FontSize',8)

text(0,0.40, ...
    '\DeltaAIC Q-L > 0 means quadratic fits better than simple linear.', ...
    'Parent',axL,'FontSize',8)

text(0,0.22, ...
    '\DeltaAIC F-Q > 0 means full model fits better than quadratic.', ...
    'Parent',axL,'FontSize',8)

%% Export figure
exportgraphics(f,fullfile(plotfolder,'brain_behav_global_quadratic_blink.pdf'),'ContentType','vector')
exportgraphics(f,fullfile(plotfolder,'brain_behav_global_quadratic_blink.svg'),'ContentType','vector')

%% Show and save table
disp(results_table)

writetable(results_table, fullfile(plotfolder,'brain_behav_model_comparison_blink.csv'))
save(fullfile(plotfolder,'brain_behav_model_comparison_blink.mat'),'results_table')