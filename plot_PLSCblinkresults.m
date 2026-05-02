%% plot integrated TFR and topo

statdims = size(stat_mse_2group_behav.stat);
stat_mse_2group_behav.stat = reshape( ...
    stat_mse_2group_behav.results.boot_result.compare_u(:,LVsel), statdims); % bootstrap ratios

close all
f = figure;

load colormap_jetlightgray.mat
% cmap = make_YAOA_colormap('match_existing', 256, true, 0);  % match_existing lab_uniform
cmap = cbrewer('div','RdBu',256);
cmap = flipud(cmap);
% % 
% figure('Color','w');
% imagesc(linspace(-1,1,256))   % fake data
% colormap(cmap)
% axis off
% % 
% title('YA (red) ← 0 (gray) → OA (blue)')

t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

%% ===== Row 1: TFR + TOPO =====
cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.clus2plot = 1;
cfg.clussign = 'pos';
cfg.integratetype = 'trapz';
cfg.subplotsize = [2 4];
cfg.subplotind = [1 2];
cfg.parameter = 'stat';
cfg.colormap = cmap;
cfg.ylabel = 'Time scale (ms)';
cfg.titleTFR = 'Feature reliability';

ft_clusterplot3D(cfg, stat_mse_2group_behav)

% Grab row-1 axes explicitly
ax_all = findall(f,'type','axes');
ax_all = flipud(ax_all); % more intuitive order after creation
axTFR  = ax_all(1);
axTopo = ax_all(2);

h = findobj(gca, 'Type', 'line');
set(h, 'LineWidth', 0.5)

axTFR.TickDir = "out";
uistack(findobj(axTFR,'Type','line'),'top')
% TFR limits and text
if LVsel == 1
    axTFR.CLim = [-40 40];
    tfr_labs   = {'-40','0','40'};
    topo_labs  = {'-200','0','200'};
    topo_ticks = [-200 0 200];
elseif LVsel == 2
    axTFR.CLim = [-10 10];
    tfr_labs   = {'-8','0','8'};
    topo_labs  = {'-35','0','35'};
    topo_ticks = [-10 0 10];
elseif LVsel == 3
    axTFR.CLim = [-8 8];
    tfr_labs   = {'-8','0','8'};
    topo_labs  = {'-35','0','35'};
    topo_ticks = [-8 0 8];
end

s = stat_mse_2group_behav.results.s;
p_lv = stat_mse_2group_behav.results.perm_result.sprob(LVsel);

axes(axTFR)  
text(1.5,100, sprintf('LV%d p = %.3f\nExpl = %.2f', ...
    LVsel, p_lv, s(LVsel).^2/sum(s.^2)), 'FontSize', 8)

% Topo limits
if LVsel == 1
    axTopo.CLim = [-200 200];
elseif LVsel == 2
    axTopo.CLim = [-100 100];
elseif LVsel == 3
    axTopo.CLim = [-35 35];
end
axis(axTopo,'square')

%% Remove auto colorbars
drawnow
delete(findall(f,'Type','ColorBar'))

%% One shared colorbar between TFR and TOPO
colormap(axTopo, cmap)
cb = colorbar(axTopo);
cb.Units = 'normalized';
cb.Box = 'off';
cb.TickDirection = 'out';
cb.FontSize = 6;
cb.TickLength = 0.02;
cb.Limits = axTopo.CLim;
cb.Ticks = [cb.Limits(1) 0 cb.Limits(2)];
cb.TickLabels = {'','',''};

% Manually place between TFR and TOPO
% cb.Position = [0.47 0.705 0.012 0.14];
cb.Position = [0.57 0.6 0.012 0.1];
% clim = caxis;   % or use your known limits
% cb.Ticks = [clim(1) 0 clim(2)];

drawnow
add_dual_colorbar_labels_figure(f, cb, tfr_labs, topo_labs);

%% ===== Row 2: BAR plot must explicitly go into tile 3 spanning both cols =====
axBar = nexttile(3,[1 2]);
hold(axBar,'on')

% Example:
% R = [r_YA; r_OA]';
% b = bar(axBar, R, 'grouped', 'BarWidth', 0.75);

% If you want only YA/OA:
R = [r_YA; r_OA]';
b = bar(axBar, R, 'grouped', 'BarWidth', 0.75);

col_YA_bar = [0.9 0.45 0.45];
col_OA_bar = [0.45 0.65 0.9];
col_YA_err = [0.75 0.15 0.15];
col_OA_err = [0.15 0.35 0.75];

b(1).FaceColor = col_YA_bar; b(1).EdgeColor = 'none';
b(2).FaceColor = col_OA_bar; b(2).EdgeColor = 'none';

xYA = b(1).XEndPoints;
xOA = b(2).XEndPoints;

errorbar(axBar, xYA, r_YA, errY_low, errY_high, ...
    'Color', col_YA_err, 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 6);
errorbar(axBar, xOA, r_OA, errO_low, errO_high, ...
    'Color', col_OA_err, 'LineStyle', 'none', 'LineWidth', 1.0, 'CapSize', 6);

for i = 1:numel(r_YA)
    if sig_labels(i) ~= ""
        ysig = max([r_YA(i)+errY_high(i), r_OA(i)+errO_high(i)]) + 0.05;
        text(axBar, mean([xYA(i) xOA(i)]), ysig, sig_labels(i), ...
            'HorizontalAlignment', 'center', ...
            'FontSize', 11, ...
            'FontWeight', 'bold');
    end
end

yline(axBar, 0, 'k')
ylim(axBar, [-1 1.1])
xlim(axBar, [0.5 numel(xYA)+0.5])
xticks(axBar, 1:numel(xYA))
xticklabels(axBar, regime_labels)
ylabel(axBar, 'Pearson''s r')
legend(axBar, {'Young','Older'}, 'Location', 'northeast')
legend(axBar,'boxoff')
box(axBar,'off')
set(axBar,'FontSize',9)

%% Figure size
w = 8.5;
h = 8.8;

f.Units = 'centimeters';
f.Position(3:4) = [w h];
f.PaperUnits = 'centimeters';
f.PaperSize = [w h];
f.PaperPosition = [0 0 w h];

% exportgraphics(f, fullfile(plotfolder,'LV1_BehavPLS.pdf'), ...
%     'ContentType','vector', 'BackgroundColor','white')
% exportgraphics(f, fullfile(plotfolder,'LV1_BehavPLS.png'), ...
%     'Resolution',300, 'BackgroundColor','white')


% plot scatters
plotscatter = 0;
if plotscatter
  zscore_scores = 1;
  nexttile
  hold on

  cols = [1 0 0;          % Young = red
    0 0.447 0.741]; % Older = blue

  r = nan(1,2);
  h_scatter = gobjects(1,2);

  for igroup = 1:2

    if zscore_scores
      behav_score = zscore(stat_mse_2group_behav.behavscores{igroup}(:,LVsel));
      brain_score = zscore(stat_mse_2group_behav.brainscores{igroup}(:,LVsel));
    else
      behav_score = stat_mse_2group_behav.behavscores{igroup}(:,LVsel);
      brain_score = stat_mse_2group_behav.brainscores{igroup}(:,LVsel);
    end

    r(igroup) = corr(behav_score, brain_score, 'rows', 'complete');

    % scatter (store handle)
    h_scatter(igroup) = scatter( brain_score, behav_score, 30, 'filled', ...
      'MarkerFaceColor', cols(igroup,:), ...
      'MarkerEdgeColor', [1 1 1], ...
      'LineWidth', 0.5);

    % regression line (not added to legend)
    p = polyfit(behav_score, brain_score, 1);
    xfit = linspace(min(behav_score), max(behav_score), 100);
    yfit = polyval(p, xfit);

    plot(xfit, yfit, '-', 'Color', cols(igroup,:), 'LineWidth', 2);
  end

  ylabel('Behavior score (z)')
  xlabel('Brain score (z)')

  % legend(h_scatter, {'Young','Older'}, 'Location', 'best')
  l = legend(h_scatter, ...
    {sprintf('YA, r = %.2f', r(1)), ...
    sprintf('OA, r = %.2f', r(2))}, ...
    'Location','best');
  % l.Position = [0.7630 0.5907 0.1389 0.1017];
  % legend boxoff
  axis padded
  set(gca, 'LineWidth', 0.5)
  box off

  title(sprintf('Brain vs Behavior\nLatent Level '), 'FontSize',9)
end
