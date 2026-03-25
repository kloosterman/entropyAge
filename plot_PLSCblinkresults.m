%% plot integrated TFR and topo
close all

f=figure;

load colormap_jetlightgray.mat

t = tiledlayout(3,3);

% f.Position = [451   695   607   254];

cfg=[];
cfg.layout = 'EEG1010.lay';
cfg.clus2plot = 1;
cfg.clussign = 'pos';
cfg.integratetype = 'trapz'; % mean or trapz
cfg.subplotsize = [2 4]; % 2 rows, 2 topo/TFR couples
cfg.subplotind = [2 3];
cfg.parameter = 'stat';
cfg.colormap = cmap;
cfg.ylabel = 'Time scale (ms)';
cfg.titleTFR = 'Feature reliability';
ft_clusterplot3D(cfg, stat_mse_2group_behav)


% plot scatters
zscore_scores = 1;
nexttile(1)
hold on

cols = [1 0 0;          % Young = red
        0 0.447 0.741]; % Older = blue

r = nan(1,2);
h_scatter = gobjects(1,2);

for igroup = 1:2

  if zscore_scores
    behav_score = zscore(stat_mse_2group_behav.behavscores{igroup}(:,1));
    brain_score = zscore(stat_mse_2group_behav.brainscores{igroup}(:,1));
  else
    behav_score = stat_mse_2group_behav.behavscores{igroup}(:,1);
    brain_score = stat_mse_2group_behav.brainscores{igroup}(:,1);
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

% title(sprintf('Young: r = %.2f \n Older: r = %.2f', r(1), r(2)))
p_lv = stat_mse_2group_behav.results.perm_result.sprob(1);
title(sprintf('Latent level, p = %.3f', p_lv), 'FontSize',9)

box off
f = gcf;

f.Units = 'centimeters';
f.Position(3:4) = [13 13];   % figure size: 12 × 8 cm

f.PaperUnits = 'centimeters';
f.PaperSize = [13 8];
f.PaperPosition = [0 0 13 8];
% set(gcf, 'Renderer', 'painters')

% exportgraphics(f,fullfile(plotfolder, 'PLSC_blink.pdf'),'ContentType','vector', 'BackgroundColor', 'white')
% exportgraphics(f,fullfile(plotfolder, 'PLSC_blink.png'),'Resolution',300,  'BackgroundColor', 'white')





% %% plot bar plots Young and Older brain_score Score vs Behavior
% 
% behav_YA = stat_mse_2group_behav.results.stacked_behavdata(1:20,:);
% behav_OA = stat_mse_2group_behav.results.stacked_behavdata(21:end,:);
% 
% brainscore_YA = stat_mse_2group_behav.brainscores{1}(:,1);
% brainscore_OA = stat_mse_2group_behav.brainscores{2}(:,1);
% 


