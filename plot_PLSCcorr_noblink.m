clear all
folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

load stat_mse_2group_noblink.mat

load colormap_jetlightgray.mat
set(groot, 'defaultFigureColor', 'w')
set(groot, 'defaultAxesColor', 'w')
set(groot, 'defaultAxesXColor', 'k')
set(groot, 'defaultAxesYColor', 'k')
set(groot, 'defaultTextInterpreter', 'none')
set(groot, 'defaultAxesTickLabelInterpreter', 'none')
%% Brain–behavior correlations: Young vs Older

% YA rows = 1:20
% OA rows = 21:end
% CI rows: 1:11 = YA, 12:22 = OA

%% Data
behav_YA = stat_mse_2group.results.stacked_behavdata(1:20,:);
behav_OA = stat_mse_2group.results.stacked_behavdata(21:end,:);

brainscore_YA = stat_mse_2group.brainscores{1}(:,1);
brainscore_OA = stat_mse_2group.brainscores{2}(:,1);

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

%% Reorder according to PLS structure
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
%% Colors (desaturated bars so difference points pop)
%% Colors
col_YA_bar   = [0.9 0.45 0.45];   % light red
col_OA_bar   = [0.45 0.65 0.9];   % light blue

col_YA_err   = [0.75 0.15 0.15];  % darker red
col_OA_err   = [0.15 0.35 0.75];  % darker blue

col_diff     = [0.15 0.15 0.15];

%% Plot
close all

R = [r_YA; r_OA]';

figure
set(gcf,'Units','centimeters')
set(gcf,'Position',[5 5 13 7.5])
set(gcf, 'InvertHardcopy', 'off')

hold on

b = bar(R,'grouped','BarWidth',0.65);

b(1).FaceColor = col_YA_bar;
b(1).EdgeColor = 'none';

b(2).FaceColor = col_OA_bar;
b(2).EdgeColor = 'none';

% b(1).EdgeColor = col_YA_err;
% b(2).EdgeColor = col_OA_err;
% b(1).LineWidth = 0.5;
% b(2).LineWidth = 0.5;

xYA = b(1).XEndPoints;
xOA = b(2).XEndPoints;

xDiff = (xYA+xOA)/2;

%% Errorbars
errorbar(xYA,r_YA,errY_low,errY_high,...
    'Color',col_YA_err,...
    'LineStyle','none',...
    'LineWidth',1.0,...
    'CapSize',6)

errorbar(xOA,r_OA,errO_low,errO_high,...
    'Color',col_OA_err,...
    'LineStyle','none',...
    'LineWidth',1.0,...
    'CapSize',6)

%% Difference points
errorbar(xDiff,r_diff,errD_low,errD_high,...
    'o',...
    'Color',col_diff,...
    'MarkerFaceColor',col_diff,...
    'MarkerEdgeColor','w',...
    'MarkerSize',7,...
    'LineWidth',1,...
    'CapSize',6)

%% Domain separators
yl = [-1 1.75];

sep = [5.5 6.5 9.5 10.5];

for s = sep
    plot([s s],yl,'Color',[.85 .85 .85],'LineWidth',1)
end

%% Significance markers
for i=1:nBehav
    if sig_labels(i)~=""
        text(xDiff(i),ul_diff(i)+.08,sig_labels(i),...
            'HorizontalAlignment','center',...
            'FontSize',11,'FontWeight','bold')
    end
end

%% Domain titles
y_domain = yl(2) + 0.1;

text(3,   y_domain,'Verbal memory','HorizontalAlignment','center','FontSize',9)
text(6,   y_domain,'Attention','HorizontalAlignment','center','FontSize',9)
text(8,   y_domain,'Working memory','HorizontalAlignment','center','FontSize',9)
text(10,  y_domain,'Fl. IQ','HorizontalAlignment','center','FontSize',9)
text(11,  y_domain,'Cr. IQ','HorizontalAlignment','center','FontSize',9)

%% Axes
yline(0,'k')

ylim(yl)
xlim([0.5 nBehav+.5])

xticks(1:nBehav)
xticklabels(pretty_labels)
xtickangle(40)

ylabel('Correlation (r) and ∆r')

legend({'Young','Older'},'Location','Best')
legend boxoff

box off
set(gca,'FontSize',9)

%% Export
% exportgraphics(gcf,'brain_behavior_age_differences.pdf','ContentType','vector')
ax = gca;
ax.Color = 'w';
ax.XColor = 'k';
ax.YColor = 'k';

set(gcf,'Units','centimeters')
set(gcf,'Position',[5 5 15 7])   % [x y width height]

% exportgraphics(gcf, 'figure.pdf', 'BackgroundColor', 'white')
exportgraphics(gcf,fullfile(plotfolder, 'brain_behav_agediff_noblink.pdf'),'ContentType','vector', 'BackgroundColor', 'white')
% print(gcf, fullfile(plotfolder, 'brain_behav_agediff_noblink.pdf'), '-dpdf', '-painters')