
%% plot 
close all

load colormap_jetlightgray.mat

fig=figure;
tiledlayout(2,3)
fig.Position = [451   695   607   254];

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

ft_clusterplot3D(cfg, stat_mse_2group_task)

%% export and save
fig.Units = 'centimeters';
fig.Position(3:4) = [13 8];   % figure size: 12 × 8 cm

fig.PaperUnits = 'centimeters';
fig.PaperSize = [13 8];
fig.PaperPosition = [0 0 13 8];

exportgraphics(fig,fullfile(plotfolder, 'taskPLS_YAvsOA.pdf'),'ContentType','vector')

%% plot mMSE time courses
close all
cfg=[];
cfg.frequency = [40 100];
cfg.avgoverfreq = 'yes';
young_mse_all_slow = ft_selectdata(cfg, young_mse_all);
old_mse_all_slow = ft_selectdata(cfg, old_mse_all);

cfg=[];
cfg.frequency = 20;
young_mse_all_fast = ft_selectdata(cfg, young_mse_all);
old_mse_all_fast = ft_selectdata(cfg, old_mse_all);
old_mse_all_fast.freq = 70;
young_mse_all_fast.freq = 70;

cfg=[];
cfg.layout = 'EEG1010.lay';
cfg.interactive = 'yes';
cfg.xlim = [-1 1];
% cfg.ylim = [1 1.25];
ft_multiplotER(cfg, old_mse_all_slow, old_mse_all_fast, young_mse_all_slow, young_mse_all_fast)
% ft_multiplotER(cfg, old_mse_all_fast, young_mse_all_fast)

%% plot posterior 20 ms, anterior 40-100 ms mse by hand
anterior = {'Fp1','Fp2','AF7','AF3','AFz','AF4','AF8', ...
            'F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
            'FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8', ...
            'T7','C5','C3','C1','Cz','C2','C4','C6','T8'};

posterior = {'TP7','CP5','CP3','CP1','CPz','CP2','CP4','CP6','TP8', ...
             'P7','P5','P3','P1','Pz','P2','P4','P6','P8', ...
             'PO7','PO3','POz','PO4','PO8','O1','Oz','O2'};

cfg = [];
cfg.channel = posterior;
cfg.avgoverchan = 'yes';
cfg.frequency = 20;
old_fast_post = ft_selectdata(cfg, old_mse_all);

cfg = [];
cfg.channel = anterior;
cfg.avgoverchan = 'yes';
cfg.frequency = [40 100];
cfg.avgoverfreq = 'yes';
old_slow_ant = ft_selectdata(cfg, old_mse_all);

cfg = [];
cfg.channel = posterior;
cfg.avgoverchan = 'yes';
cfg.frequency = 20;
young_fast_post = ft_selectdata(cfg, young_mse_all);

cfg = [];
cfg.channel = anterior;
cfg.avgoverchan = 'yes';
cfg.frequency = [40 100];
cfg.avgoverfreq = 'yes';
young_slow_ant = ft_selectdata(cfg, young_mse_all);

% subject x time
dat_old_fast_post   = squeeze(old_fast_post.powspctrm);
dat_old_slow_ant    = squeeze(old_slow_ant.powspctrm);
dat_young_fast_post = squeeze(young_fast_post.powspctrm);
dat_young_slow_ant  = squeeze(young_slow_ant.powspctrm);

time = young_slow_ant.time;

% means
m_old_fast_post   = mean(dat_old_fast_post, 1);
m_old_slow_ant    = mean(dat_old_slow_ant, 1);
m_young_fast_post = mean(dat_young_fast_post, 1);
m_young_slow_ant  = mean(dat_young_slow_ant, 1);

% SEM
sem_old_fast_post   = std(dat_old_fast_post, 0, 1) ./ sqrt(size(dat_old_fast_post, 1));
sem_old_slow_ant    = std(dat_old_slow_ant, 0, 1) ./ sqrt(size(dat_old_slow_ant, 1));
sem_young_fast_post = std(dat_young_fast_post, 0, 1) ./ sqrt(size(dat_young_fast_post, 1));
sem_young_slow_ant  = std(dat_young_slow_ant, 0, 1) ./ sqrt(size(dat_young_slow_ant, 1));

f=figure;
f.Position = [    1     1   438   282]
     
hold on

% Colors (Old = blue, Young = red)
c_old   = [0.0 0.3 0.8];
c_young = [0.8 0.1 0.1];

% Shaded error bands
fill([time fliplr(time)], ...
     [m_old_fast_post + sem_old_fast_post, fliplr(m_old_fast_post - sem_old_fast_post)], ...
     c_old, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

fill([time fliplr(time)], ...
     [m_old_slow_ant + sem_old_slow_ant, fliplr(m_old_slow_ant - sem_old_slow_ant)], ...
     c_old, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

fill([time fliplr(time)], ...
     [m_young_fast_post + sem_young_fast_post, fliplr(m_young_fast_post - sem_young_fast_post)], ...
     c_young, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

fill([time fliplr(time)], ...
     [m_young_slow_ant + sem_young_slow_ant, fliplr(m_young_slow_ant - sem_young_slow_ant)], ...
     c_young, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');

% Mean lines
h1 = plot(time, m_old_fast_post,   '--', 'Color', c_old,   'LineWidth', 2); % dashed
h2 = plot(time, m_old_slow_ant,    '-',  'Color', c_old,   'LineWidth', 2);

h3 = plot(time, m_young_fast_post, '--', 'Color', c_young, 'LineWidth', 2); % dashed
h4 = plot(time, m_young_slow_ant,  '-',  'Color', c_young, 'LineWidth', 2);

xlim([-1 1])
xline(0)
xlabel('Time (s)')
ylabel('Sample entropy')

l = legend([h1 h2 h3 h4], ...
       {'Older fast post.','Older slow ant.', ...
        'Young fast post.','Young slow ant.'}, ...
       'Location','best');
l.Position = [   0.5708    0.7980    0.3311    0.2181];
box off
%% Export
% exportgraphics(gcf,'brain_behavior_age_differences.pdf','ContentType','vector')
set(gcf,'Units','centimeters')
set(gcf,'Position',[5 5 10 7])   % [x y width height]

exportgraphics(gcf,fullfile(plotfolder, 'MSEtimecourses.pdf'),'ContentType','vector')
