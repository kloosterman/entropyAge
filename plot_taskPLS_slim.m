%% Extract brain scores
Y = stat_mse_2group_task.brainscores{1}(:,1);
O = stat_mse_2group_task.brainscores{2}(:,1);

% Combine + z-score across all subjects
B = [Y; O];
Bz = zscore(B);

% Split back
nY = numel(Y);
Y = Bz(1:nY);
O = Bz(nY+1:end);

%% Remove outliers |z| > 3
outY = abs(Y) > 3;
outO = abs(O) > 3;

fprintf('Brain score outliers removed: Young=%d, Older=%d\n', sum(outY), sum(outO));

Y(outY) = NaN;
O(outO) = NaN;

%% Stats
meanY = mean(Y, 'omitnan');
meanO = mean(O, 'omitnan');

semY = std(Y, 'omitnan') / sqrt(sum(~isnan(Y)));
semO = std(O, 'omitnan') / sqrt(sum(~isnan(O)));

[~, p] = ttest2(Y, O);

%% Plot
% figure('Units','centimeters','Position',[5 5 5 6]);
nexttile(4)
hold on

col_Y = [1 0 0];
col_O = [0 0.447 0.741];
jitter = 0.08;

% Bars
x = 1.5; % single category

b = bar(x, [meanY meanO], 0.8); % 0.8 = total group width
hold on

b(1).FaceColor = col_Y;
b(2).FaceColor = col_O;

b(1).FaceAlpha = 0.6;
b(2).FaceAlpha = 0.6;
set(b, 'EdgeColor', 'none')

% Scatter
% Get bar centers
xY = b(1).XEndPoints;
xO = b(2).XEndPoints;

% Jitter around each bar
x_jit_Y = xY + (rand(sum(~isnan(Y)),1)-0.5)*2*jitter;
x_jit_O = xO + (rand(sum(~isnan(O)),1)-0.5)*2*jitter;

scatter(x_jit_Y, Y(~isnan(Y)), 4, col_Y, 'filled', ...
    'MarkerFaceAlpha',0.4, 'MarkerEdgeColor','w','LineWidth',0.25);

scatter(x_jit_O, O(~isnan(O)), 4, col_O, 'filled', ...
    'MarkerFaceAlpha',0.4, 'MarkerEdgeColor','w','LineWidth',0.25);

% Error bars
errorbar(xY, meanY, semY, 'k', 'LineStyle','none','LineWidth',0.5);
errorbar(xO, meanO, semO, 'k', 'LineStyle','none','LineWidth',0.5);

%% Significance
ymax = max([meanY+semY, meanO+semO]);
ymin = min([meanY-semY, meanO-semO]);
yr = ymax - ymin;
if yr == 0, yr = 1; end

% y_sig = ymax + 0.15*yr;
y_sig = 2.25;

if p < 0.001
    stars = '***';
elseif p < 0.01
    stars = '**';
elseif p < 0.05
    stars = '*';
else
    stars = 'n.s.';
end

plot([1.25 1.75], [y_sig y_sig], 'k', 'LineWidth',1);
text(1.5, y_sig + 0.05*yr, stars, 'HorizontalAlignment','center');

%% Formatting
xticks(1.5)
% xticklabels({'blink-aligned rs-EEG entropy'})
% xtickangle(25)

% Remove x tick labels
% xticks([])
xticklabels({})
% Add manual centered label
text(1.5, ax.YLim(1) - 0.025*range(ax.YLim), ...
    {'Brain score'}, ... ,'Young vs. Older'
    'HorizontalAlignment','right', ...
    'VerticalAlignment','top', ...
    'FontSize',8, ...
    'Rotation',25, ...
    'Clipping','off');
box off

ax = gca;
ax.FontSize = 8;
ax.TickDir = 'out';
ax.XLim = [1.2 1.8];
if isempty(oa_excl)
  ax.YLim = [-3.4 3];
else
  ax.YLim = [-3 3];
end
ax.YAxis.Visible = 'off';
ax.XAxis.FontSize = 7;
pos = ax.Position;
pos(3) = pos(3) * 1.3;   % increase width (try 1.2–1.5)
ax.Position = pos;

lgd = legend({'Young','Older'}, 'Location', 'best');
lgd.ItemTokenSize = [5 5];
lgd.Position = [0.4825    0.7535    0.1693    0.0580];
lgd.TextColor = 'k';
legend boxoff

%% topo
stat = stat_mse_2group_task;

% Select scales: 40–100 ms
freqsel = stat.freq >= 20 & stat.freq <= 100;

% Average stat over selected scales and all time bins
topo_avg = squeeze(mean(mean(stat.stat(:, freqsel, :), 3, 'omitnan'), 2, 'omitnan'));
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
mask_avg = squeeze(mean(mean(stat.posclusterslabelmat(:, freqsel, :), 3), 2)) > 0;

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
cfg.gridscale = 300;   % default ~67 → increase a lot
t = title('Brain signal variability', ...
    'HorizontalAlignment','left');
% t.Position(1) = 0.4;
% t.Position(2) = 0.85;
t2 = subtitle('Young vs. Older', ...
    'HorizontalAlignment','left');

nexttile(5,[1 2]); 
ft_topoplotER(cfg, topo);
cmap = cbrewer('div','RdBu',256);
cmap = flipud(cmap);
colormap(cmap)
axis tight
h = findall(gca, 'Type', 'line');
set(h, 'LineWidth', 0.5);
% title('Task PLS stat, 40–100 ms scales, all time bins');
% t = title('Brain signal variability', ...
%     'HorizontalAlignment','right');
% t2 = subtitle('Young vs. Older', ...
%     'HorizontalAlignment','right');
% t.Position(1) = 0.4;
% t.Position(2) = 0.85;
% h=annotation('textbox', [ 0.7031    0.9240    0.2899    0.0432], 'String', 'Young vs. Older', ...
%            'EdgeColor','none','FontSize',8);
cb = colorbar;
cb.Location = 'south';   % horizontal, anchored at bottom of tile
cb.Position =   [ 0.82    0.75    0.1    0.02];
cb.FontSize = 7;
cb.Box = 'off';
xlabel(cb, 'Bootstrap ratio')

ax = gca;
clim = ax.CLim;