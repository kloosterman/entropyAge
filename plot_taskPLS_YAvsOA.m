
%% plot 
close all

% load colormap_jetlightgray.mat
cmap = cbrewer('div','RdBu',256);
cmap = flipud(cmap);

fig=figure;
tiledlayout(2,4,'TileSpacing','compact','Padding','compact')
fig.Position = [451   695   800   254];

cfg = [];
cfg.layout = 'EEG1010.lay';
cfg.clus2plot = 1;
cfg.clussign = 'pos';
cfg.integratetype = 'mean';
cfg.parameter = 'stat';
cfg.colormap = cmap;
cfg.subplotind = [1 2 3];
cfg.topo_usemask = 'no';
% cfg.TFRspan = [1 2];        % TFR spans 2 tiles

cfg.topo_windows(1).timerange = [-0.8 0.4];
cfg.topo_windows(1).scalerange = [0 20];
cfg.topo_windows(1).title = 'Fast scales';

cfg.topo_windows(2).timerange = [-0.8 0.4];
cfg.topo_windows(2).scalerange = [60 100];
cfg.topo_windows(2).title = 'Slow scales';

ft_clusterplot3D(cfg, stat_mse_2group_task)
nexttile(2); axis tight; set(findall(gcf,'-property','LineWidth'),'LineWidth',0.5)
nexttile(3); axis tight; set(findall(gcf,'-property','LineWidth'),'LineWidth',0.5)


%%
stat = stat_mse_2group_task;

cfg_tc = [];
cfg_tc.layout = 'EEG1010.lay';
cfg_tc.mask = stat.posclusterslabelmat;                    % chan x freq x time
cfg_tc.freqrange = [20 20];
cfg_tc.timerange = [-1.25 1.25];
cfg_tc.weighting = 'cluster';
cfg_tc.splitmode = 'medianY';

out20 = extract_ant_post_timecourses(cfg_tc, stat, young_mse_all.powspctrm, old_mse_all.powspctrm);
cfg_tc.freqrange = [40 100];
out60 = extract_ant_post_timecourses(cfg_tc, stat, young_mse_all.powspctrm, old_mse_all.powspctrm);

% figure; hold on
nexttile(4); hold on; axis tight;

% Colors
% col_YA = [1 0 0]; % red
% col_OA = [0 0 1]; % blue
col_YA = [0.85 0.2 0.2];   % slightly softer red
col_OA = [0.2 0.4 0.8];    % nicer blue (less saturated than [0 0 1])

% --- Compute means and SEM ---
% YA anterior (60)
m_YA_ant = mean(out60.ant.YA,1);
sem_YA_ant = std(out60.ant.YA,[],1) ./ sqrt(size(out60.ant.YA,1));

% YA posterior (20)
m_YA_post = mean(out20.post.YA,1);
sem_YA_post = std(out20.post.YA,[],1) ./ sqrt(size(out20.post.YA,1));

% OA anterior (60)
m_OA_ant = mean(out60.ant.OA,1);
sem_OA_ant = std(out60.ant.OA,[],1) ./ sqrt(size(out60.ant.OA,1));

% OA posterior (20)
m_OA_post = mean(out20.post.OA,1);
sem_OA_post = std(out20.post.OA,[],1) ./ sqrt(size(out20.post.OA,1));

t = stat.time;

% --- Helper for shaded error ---
plot_shaded = @(x, m, sem, col, ls) ...
    fill([x fliplr(x)], [m-sem fliplr(m+sem)], col, ...
    'FaceAlpha', 0.2, 'EdgeColor', 'none');

% --- Plot shaded areas ---
% --- Plot shaded areas (HIDE from legend) ---
h1 = plot_shaded(t, m_YA_ant, sem_YA_ant, col_YA);
h2 = plot_shaded(t, m_YA_post, sem_YA_post, col_YA);
h3 = plot_shaded(t, m_OA_ant, sem_OA_ant, col_OA);
h4 = plot_shaded(t, m_OA_post, sem_OA_post, col_OA);

set([h1 h2 h3 h4], 'HandleVisibility', 'off');
% --- Plot lines on top ---
hYA_ant = plot(t, m_YA_ant, '-',  'Color', col_YA, 'LineWidth', 1.5);
hYA_post = plot(t, m_YA_post, ':', 'Color', col_YA, 'LineWidth', 1.5);

hOA_ant = plot(t, m_OA_ant, '-',  'Color', col_OA, 'LineWidth', 1.5);
hOA_post = plot(t, m_OA_post, ':', 'Color', col_OA, 'LineWidth', 1.5);

% --- Decorations ---
xline(0,'k:')
ax2 = nexttile(8);
axis off  % hide axes
lgd = legend(ax2, [hYA_ant hYA_post hOA_ant hOA_post], ...
       {'YA ant60','YA post20','OA ant60','OA post20'}, ...
       'Location','southoutside', ...
       'Orientation','vertical');
lgd.Position =  [  0.7603    0.3047    0.1974    0.1595];
legend boxoff
xlabel('Time (s)')
ylabel('Entropy')
box on

%% Export
% exportgraphics(gcf,'brain_behavior_age_differences.pdf','ContentType','vector')
set(gcf,'Units','centimeters')
set(gcf,'Position',[5 5 13 7])   % [x y width height]

% export and save
fig.Units = 'centimeters';
fig.Position(3:4) = [13 7];   % figure size: 12 × 8 cm

fig.PaperUnits = 'centimeters';
fig.PaperSize = [13 7];
fig.PaperPosition = [0 0 13 7];

% exportgraphics(fig,fullfile(plotfolder, 'taskPLS_YAvsOA.pdf'),'ContentType','vector')

% exportgraphics(gcf,fullfile(plotfolder, 'taskPLS_corrbehavior.pdf'),'ContentType','vector')
% exportgraphics(gcf,fullfile(plotfolder, 'taskPLS_corrbehavior.png'),'ContentType','vector')

