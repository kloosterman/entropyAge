clear
folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

load stat_mse_2group_noblink.mat

load colormap_jetlightgray.mat

%% plot integrated topo
close all

f=figure;
tiledlayout(2,3)
f.Position = [451   695   607   254];

cfg=[];
cfg.layout = 'EEG1010.lay';
cfg.clus2plot = 1;
cfg.clussign = 'pos';
cfg.integratetype = 'trapz'; % mean or trapz
cfg.subplotsize = [2 4]; % 2 rows, 2 topo/TFR couples
cfg.subplotind = [1 2];
cfg.parameter = 'stat';
cfg.colormap = cmap;
cfg.ylabel = 'Time scale (ms)';
ft_clusterplot3D(cfg, stat_mse_2group)

nexttile(1);
xlim([15 105])
ylim([-5 50])
set(gca, 'XTick', 0:20:100);
set(gca, 'YTick', 0:20:40);
box off

%% plot scatters

% f = figure; 
% f.Position = [100 100 500 400];
nexttile
hold on

cols = [1 0 0;          % Young = red
        0 0.447 0.741]; % Older = blue

r = nan(1,2);
h_scatter = gobjects(1,2);

for igroup = 1:2

    behav = zscore(stat_mse_2group.behavscores{igroup}(:,1));
    brain = zscore(stat_mse_2group.brainscores{igroup}(:,1));

    r(igroup) = corr(behav, brain, 'rows', 'complete');

    % scatter (store handle)
    h_scatter(igroup) = scatter(behav, brain, 30, 'filled', ...
        'MarkerFaceColor', cols(igroup,:), ...
        'MarkerEdgeColor', [1 1 1], ...
        'LineWidth', 0.5);

    % regression line (not added to legend)
    p = polyfit(behav, brain, 1);
    xfit = linspace(min(behav), max(behav), 100);
    yfit = polyval(p, xfit);

    plot(xfit, yfit, '-', 'Color', cols(igroup,:), 'LineWidth', 2);
end

xlabel('Behavior score (z)')
ylabel('Brain score (z)')

% legend(h_scatter, {'Young','Older'}, 'Location', 'best')
l = legend(h_scatter, ...
    {sprintf('YA, r = %.2f', r(1)), ...
     sprintf('OA, r = %.2f', r(2))}, ...
    'Location','best');
l.Position = [0.7630 0.5907 0.1389 0.1017];
% legend boxoff
axis padded

% title(sprintf('Young: r = %.2f \n Older: r = %.2f', r(1), r(2)))

box off
set(gca,'LineWidth',1)

fig = gcf;

fig.Units = 'centimeters';
fig.Position(3:4) = [13 8];   % figure size: 12 × 8 cm

fig.PaperUnits = 'centimeters';
fig.PaperSize = [13 8];
fig.PaperPosition = [0 0 13 8];

exportgraphics(fig,fullfile(plotfolder, 'PLSC_noblink.pdf'),'ContentType','vector', 'BackgroundColor', 'white')

% %% plot bar plots Young and Older Brain Score vs Behavior
% 
% behav_YA = stat_mse_2group.results.stacked_behavdata(1:20,:);
% behav_OA = stat_mse_2group.results.stacked_behavdata(21:end,:);
% 
% brainscore_YA = stat_mse_2group.brainscores{1}(:,1);
% brainscore_OA = stat_mse_2group.brainscores{2}(:,1);
% 


