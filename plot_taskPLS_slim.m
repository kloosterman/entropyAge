%% Plot task PLS bootstrap-ratio topo and time-by-entropy map
figLayout = findobj(gcf, 'Type', 'tiledlayout');
figLayout = figLayout(1);

%% Bootstrap-ratio volume
stat = stat_mse_2group_task;
statdims = size(stat.stat);
if isfield(stat, 'results') && isfield(stat.results, 'boot_result') && ...
        isfield(stat.results.boot_result, 'compare_u')
    bsr = reshape(stat.results.boot_result.compare_u(:,1), statdims);
else
    bsr = stat.stat;
end
stat.stat = bsr;
if isfield(stat, 'posclusterslabelmat')
    bsr_mask = logical(stat.posclusterslabelmat);
else
    bsr_mask = abs(bsr) > 3;
end

% Select scales: 20-100 ms
freqsel = stat.freq >= 20 & stat.freq <= 100;

% Average bootstrap ratios over selected scales and all time bins.
topo_avg = squeeze(mean(mean(bsr(:, freqsel, :), 3, 'omitnan'), 2, 'omitnan'));
% result: [60 x 1]

% Put into FieldTrip-like structure for topoplot
topo = [];
topo.label    = stat.label;
topo.avg      = topo_avg;
topo.dimord   = 'chan';
topo.time     = 0;   % dummy
topo.label    = stat.label;

% % Plot
% Average positive cluster mask over same dimensions
mask_avg = squeeze(mean(mean(bsr_mask(:, freqsel, :), 3), 2)) > 0;

topo.mask = mask_avg;

cfg = [];
cfg.layout          = 'EEG1010.lay';
cfg.parameter       = 'avg';
% cfg.maskparameter   = 'mask';
cfg.maskstyle       = 'outline';
cfg.colorbar        = 'no';
cfg.comment         = 'no';
cfg.marker          = 'off';
cfg.zlim = 'maxabs';
cfg.figure = 'gca';
cfg.gridscale = 300;   % default ~67 -> increase a lot

tileAx = nexttile(figLayout, 4, [1 3]);
tilePos = tileAx.Position;
tileAx.Visible = 'off';
tileAx.XTick = [];
tileAx.YTick = [];
tileAx.Color = 'none';
box(tileAx, 'off')

gap = 0.014;
tfrPos = [tilePos(1) + tilePos(3)*0.12, ...
          tilePos(2) + tilePos(4)*0.06, ...
          tilePos(3)*0.58, ...
          tilePos(4)*0.42];
topoPos = [tilePos(1) + tilePos(3)*0.18, ...
           tilePos(2) + tilePos(4)*0.52 + gap, ...
           tilePos(3)*0.56, ...
           tilePos(4)*0.46];
clim_bsr = [-3 3];

topoAx = axes('Position', topoPos);
ft_topoplotER(cfg, topo);
cmap = cbrewer('div','RdBu',256);
cmap = flipud(cmap);
colormap(cmap)
caxis(clim_bsr)
axis tight
h = findall(gca, 'Type', 'line');
set(h, 'LineWidth', 0.5);
title('')
cb = colorbar;
cb.Location = 'eastoutside';
cb.Position = [topoPos(1) + topoPos(3) + 0.006, ...
               topoPos(2) + topoPos(4)*0.20, ...
               0.008, ...
               topoPos(4)*0.60];
cb.FontSize = 7;
cb.Box = 'off';
cb.Label.String = '';

%% Time-by-entropy BSR map
tfrAx = axes('Position', tfrPos);
tfr_bsr = squeeze(mean(bsr, 1, 'omitnan'));
imagesc(tfrAx, stat.time, stat.freq, tfr_bsr);
set(tfrAx, 'YDir', 'normal')
colormap(tfrAx, cmap)
caxis(tfrAx, clim_bsr)
hold(tfrAx, 'on')
xline(tfrAx, 0, 'k:', 'LineWidth', 0.5, 'HandleVisibility', 'off')

xlabel(tfrAx, 'Time (s)', 'FontSize', 7)
ylabel(tfrAx, 'Time scale (ms)', 'FontSize', 6)
box(tfrAx, 'off')
tfrAx.FontSize = 7;
tfrAx.TickDir = 'out';
dt = median(diff(stat.time));
df = median(diff(stat.freq));
tfrAx.XLim = [min(stat.time)-dt/2 max(stat.time)+dt/2];
tfrAx.YLim = [min(stat.freq)-df/2 max(stat.freq)+df/2];
tfrAx.XTick = [-1 0 1];
tfrAx.YTick = stat.freq;
