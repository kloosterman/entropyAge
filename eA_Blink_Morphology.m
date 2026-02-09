%% eA_Blink_feature_extraction
% Refactored script for extracting blink features for older and younger adults

%% default settings
restoredefaultpath;
addpath(fullfile('C:','Users','morit','Desktop','Toolboxes_MATLAB','fieldtrip-20240916'));
ft_defaults;

%% paths and generic settings
data_path = fullfile('C:','Users','morit','Desktop','FoPra_Daten','blink_morphology');
plot_individual = true;                % turn on/off per-subject figures
time = linspace(-1.5, 1.5, 1051);

%% subject lists (unchanged)
old_ids = {'8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45', '5HO51', ...
           '7MT51', '4CH63', '7CU61', '7WE49', '6RE54', '4KS63', '6JH59', ...
           '7HI40', '6CR44', '6CM48', '10SH66', '11MS53'};

young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89', '4ML96', '5SR93', '11VZ96', '10SH95'};

%% Process groups (keeps outputs separate, no overwritten variables)
old_out   = process_group(old_ids,  'OA_', data_path, time, plot_individual);
young_out = process_group(young_ids, 'YA_', data_path, time, plot_individual);

%% Example: quick comparison plot (grand averages overlayed)
figure('Color','w'); hold on;
plot(time, old_out.VEOG_gavg, 'LineWidth', 2, 'Color', [0.85 0.33 0.10]); % orange/red
plot(time, young_out.VEOG_gavg, 'LineWidth', 2, 'Color', [0 0.45 0.74]);   % blue
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
legend({'Old - grand avg','Young - grand avg'}, 'Location','Best');
title('Grand Average Blink by Group'); grid on;

%% Save results optionally
% save('blink_features_summary.mat','old_out','young_out');

%% ---------------------- Local functions -------------------------------
function out = process_group(ids, prefix, data_path, time, plot_individual)
% process_group - loads subject blink files, computes per-subject VEOG average
% and blink measures. Returns a struct with fields:
%   .IDs, .VEOG_all, .VEOG_gavg, .measures_table

nSub = numel(ids);
VEOG_all = NaN(nSub, numel(time));
M1 = NaN(nSub,1); M2 = NaN(nSub,1); M3 = NaN(nSub,1); M4 = NaN(nSub,1); M5 = NaN(nSub,1);

for i = 1:nSub
    sid = ids{i};
    filename = fullfile(data_path, [prefix sid '_blink.mat']);

    if ~exist(filename,'file')
        warning('File not found: %s', filename);
        continue;
    end

    tmp = load(filename);
    vars = fieldnames(tmp);
    data = tmp.(vars{1});   % first variable in the MAT-file

    % basic validation
    if isempty(data)
        warning('Subject %s - data empty, skipping.', sid); continue;
    end
    if size(data,2) ~= numel(time)
        warning('Subject %s - unexpected sample length (%d). Skipping.', sid, size(data,2));
        continue;
    end

    % select VEOG rows (assumed to be rows 1 and 2) and detect trials with NaNs
    data_sel = data(1:2,:,:);  % 2 x T x N
    trial_has_nan = squeeze( any( any( isnan(data_sel), 1 ), 2 ) );
    valid_idx = find(~trial_has_nan);

    if isempty(valid_idx)
        warning('Subject %s - no valid blink trials (all contain NaN). Skipping.', sid);
        continue;
    end

    % average across valid trials
    data_clean = data(:,:,valid_idx);      % 4 x T x n_valid
    avg_data = mean(data_clean,3);         % 4 x T

    % rereference VEOG and compute mean VEOG
    VEOG1 = avg_data(1,:);
    VEOG2 = -avg_data(2,:);
    VEOG_avg = (VEOG1 + VEOG2) / 2;
    VEOG_all(i,:) = VEOG_avg;

    % optional: individual plotting
    if plot_individual
        figure('Color','w'); hold on;
        plot(time, VEOG1, 'Color', [0.8 0.8 0.8]);
        plot(time, VEOG2, 'Color', [0.8 0.8 0.8]);
        plot(time, VEOG_avg, 'LineWidth', 2);
        xlabel('Time (s)'); ylabel('Amplitude');
        title(sprintf('Subject %s – Average Blink (VEOG)', sid)); grid on;
    end

    % compute measures
    [m1,m2,m3,m4,m5] = computeMeasures(VEOG_avg, time);
    M1(i)=m1; M2(i)=m2; M3(i)=m3; M4(i)=m4; M5(i)=m5;
end

% assemble outputs
out.IDs = ids(:);
out.VEOG_all = VEOG_all;
out.VEOG_gavg = nanmean(VEOG_all,1);
out.measures_table = table(out.IDs, M1, M2, M3, M4, M5, ...
    'VariableNames', {'ID','M1_peak_width','M2_peak_amp','M3_peak_to_neg','M4_neg_to_zero','M5_peak_neg_amp'});

end

function [M1,M2,M3,M4,M5] = computeMeasures(y, t)
% computeMeasures - compute M1..M5 from 1D signal y at times t
% M2: peak amplitude and index
[M2, peak_idx] = max(y);

% M1: FWHM
half_height = M2 / 2;
left_idx = find(y(1:peak_idx) <= half_height, 1, 'last');
right_rel = find(y(peak_idx:end) <= half_height, 1, 'first');
if ~isempty(right_rel), right_idx = peak_idx - 1 + right_rel; else right_idx = []; end
if isempty(left_idx) || isempty(right_idx), M1 = NaN; else M1 = t(right_idx) - t(left_idx); end

% M3: peak to post-peak negativity
if peak_idx < length(y)
    [~, neg_rel_idx] = min(y(peak_idx:end));
    neg_idx = peak_idx + neg_rel_idx - 1;
    M3 = t(neg_idx) - t(peak_idx);
else
    M3 = NaN; neg_idx = [];
end

% M4: negative trough to first zero crossing after it
if ~isempty(neg_idx)
    post_neg = y(neg_idx:end);
    zero_rel = find(post_neg >= 0, 1, 'first');
    if isempty(zero_rel)
        M4 = NaN;
    else
        zero_idx = neg_idx + zero_rel - 1;
        M4 = t(zero_idx) - t(neg_idx);
    end
else
    M4 = NaN;
end

% M5: peak negativity amplitude after positive peak
if peak_idx < length(y)
    M5 = min(y(peak_idx:end));
else
    M5 = NaN;
end

end

% %% Grand average younger and older adults in one plot
% 
% figure('Color','w'); hold on;
% 
% plot(time, young_VEOG_gavg, 'LineWidth', 2, 'Color', 'r');
% plot(time, old_VEOG_gavg, 'LineWidth', 2, 'Color', 'b');
% xlabel('Time (s)')
% ylabel('Amplitude (µV)');
% title('Grand average blink for younger and older adults');
% grid on;
% legend({'Younger adults', 'Older adults'});

%% Plot averaged blink + morphology arrows for a single OLD subject (robust)
subjID = '5MA56';

% --- parameters / paths (edit data_path if needed) ---
data_path = fullfile('C:','Users','morit','Desktop','FoPra_Daten','blink_morphology');
time = linspace(-1.5, 1.5, 1051);   % ensure this matches your data
plot_individual = false;             % not used here but kept for consistency

% --- Try to find subject in old_out struct first (if available) ------------
use_struct = exist('old_out','var') && isstruct(old_out);
if use_struct
    idx = find(strcmp(old_out.IDs, subjID), 1);
else
    % fallback to check subject file on disk
    idx = [];
end

% --- load the data either from struct or file --------------------------------
if ~isempty(idx)
    % take VEOG averaged waveform from old_out if available
    % old_out.VEOG_all is subjects x time
    if size(old_out.VEOG_all,1) >= idx && any(~isnan(old_out.VEOG_all(idx,:)))
        blink_amp = old_out.VEOG_all(idx,:);
        t = time;
    else
        error('Subject %s present in old_out but VEOG_all empty for that index.', subjID);
    end
else
    % load file OA_<ID>_blink.mat directly
    fname = fullfile(data_path, ['OA_' subjID '_blink.mat']);
    if ~exist(fname,'file')
        error('Subject file not found: %s\nEither create old_out or fix data_path.', fname);
    end
    tmp = load(fname);
    vars = fieldnames(tmp);
    data = tmp.(vars{1});   % expected: 4 x Ntime x Nblink

    if isempty(data)
        error('No data in file for subject %s.', subjID);
    end
    if size(data,2) ~= numel(time)
        warning('Time axis length differs from expected. Adapting time vector to data length.');
        t = linspace(-1.5, 1.5, size(data,2));
    else
        t = time;
    end

    % keep only VEOG channels and drop trials having any NaN in channels 1:2
    data_sel = data(1:2,:,:);
    trial_has_nan = squeeze( any( any( isnan(data_sel), 1 ), 2 ) );
    valid_idx = find(~trial_has_nan);

    if isempty(valid_idx)
        error('No valid blink trials for subject %s (all trials contain NaN).', subjID);
    end

    data_clean = data(:,:,valid_idx);          % 4 x T x n_valid
    avg_data = mean(data_clean,3);            % 4 x T average over valid trials

    VEOG1 = avg_data(1,:);
    VEOG2 = -avg_data(2,:);
    blink_amp = (VEOG1 + VEOG2) / 2;          % positive-going VEOG average
end

% --- basic plotting -------------------------------------------------------
f = figure('Color','w','Units','centimeters','Position',[2.5 2.5 16.0 8.9]);
set(f, 'PaperUnits','centimeters', 'PaperPosition',[2.5 2.5 16.0 8.9],'PaperSize',[21.0 29.7]);

plot(t, blink_amp, 'b', 'LineWidth',1.5); hold on;
xlabel('Time (s)'); ylabel('Amplitude (\muV)');
title(sprintf('Averaged Blink (VEOG) - %s', subjID));
set(gca,'FontSize',10,'Color','w'); grid off; box on;

% reference lines
yline(0,'Color',[0.3 0.3 0.3],'LineStyle',':','LineWidth',1);
xline(0,'k:','LineWidth',1.5);

% --- find positive peak and the subsequent negative trough (post t>0) -----
[posPeakVal,posPeakIdx] = max(blink_amp);
posPeakTime = t(posPeakIdx);

% negative trough after positive peak (restrict search to t>0 to match your original)
pos_mask = t > 0;
if any(pos_mask)
    local_seq = blink_amp(pos_mask);
    [negPeakVal_local, relIdx] = min(local_seq);
    tmpIdx = find(pos_mask);
    negPeakIdx = tmpIdx(relIdx);
    negPeakVal = negPeakVal_local;
    negPeakTime = t(negPeakIdx);
else
    negPeakIdx = []; negPeakVal = NaN; negPeakTime = NaN;
end

% % plot peak markers
% plot(posPeakTime, posPeakVal, 'ko', 'MarkerFaceColor','k');
% if ~isempty(negPeakIdx)
%     plot(negPeakTime, negPeakVal, 'ko', 'MarkerFaceColor','k');
% end

% dotted horizontal lines for peaks (safe bounds)
xr = xlim; yr = ylim;
line([posPeakTime min(posPeakTime+0.5,xr(2))],[posPeakVal posPeakVal],'Color',[0.3 0.3 0.3],'LineStyle',':');
if ~isnan(negPeakVal)
    line([negPeakTime min(negPeakTime+0.5,xr(2))],[negPeakVal negPeakVal],'Color',[0.3 0.3 0.3],'LineStyle',':');
end

% --- Morphology measure endpoints M1 / M3 ---------------------------------
% left zero crossing before the positive peak:
left_zero_idx = find(blink_amp(1:posPeakIdx) <= 0, 1, 'last');
% right zero crossing after positive peak:
right_rel = find(blink_amp(posPeakIdx:end) <= 0, 1, 'first');
if ~isempty(right_rel)
    right_zero_idx = posPeakIdx - 1 + right_rel;
else
    right_zero_idx = []; % handle later
end

% safe fallbacks if zero crossing not found
if isempty(left_zero_idx)
    left_zero_idx = 1;
end
if isempty(right_zero_idx)
    right_zero_idx = numel(blink_amp);
end

xL = t(left_zero_idx);
xR = t(right_zero_idx);

% M3 endpoints are posPeakTime to negPeakTime (if available), else fallback
M3_x1 = posPeakTime;
if ~isempty(negPeakIdx)
    M3_x2 = negPeakTime;
else
    M3_x2 = xR;    % fallback horizontal span
end
M3_y = min(blink_amp);  % conservative placement if trough not found

% --- annotation coordinate converter (axes -> normalized figure coords) ---
ax = gca;
xlimv = get(ax,'XLim'); ylimv = get(ax,'YLim');
axpos = get(ax,'Position');     % normalized axes position in figure

to_norm = @(xd, yd) deal( ...
    axpos(1) + ((xd - xlimv(1))/(xlimv(2)-xlimv(1))) * axpos(3), ...
    axpos(2) + ((yd - ylimv(1))/(ylimv(2)-ylimv(1))) * axpos(4) ...
);

% vertical range to shift labels slightly
yrange = diff(ylimv);
shift_M1_y = -0.02*yrange;
shift_M3_y = -0.01*yrange;
shift_M2_x = 0.5;   % seconds shift for vertical M2 arrow
shift_M4_x = 0.75;  % seconds shift for vertical M4 arrow

% --- M1: horizontal width between left & right zero crossings (shifted a bit vertically) ---
[nx1, ny1] = to_norm(xL, 0 + shift_M1_y);
[nx2, ny2] = to_norm(xR, 0 + shift_M1_y);
annotation('line',[nx1 nx2],[ny1 ny2],'LineWidth',1.4);
annotation('arrow',[nx1 nx2],[ny1 ny2],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
annotation('arrow',[nx2 nx1],[ny2 ny1],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
% place M1 text in data coordinates (so it follows axis scaling)
text(mean([xL xR]), 0 + 2*shift_M1_y, 'M_1', 'HorizontalAlignment','center','FontWeight','bold');

% --- M2: vertical arrow from 0 to positive peak (shift right by shift_M2_x) ---
posM2_x_shifted = min(max(posPeakTime + shift_M2_x, xlimv(1)), xlimv(2));
[nx1, ny1] = to_norm(posM2_x_shifted, 0);
[nx2, ny2] = to_norm(posM2_x_shifted, posPeakVal);
annotation('line',[nx1 nx2],[ny1 ny2],'LineWidth',1.4);
annotation('arrow',[nx1 nx2],[ny1 ny2],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
annotation('arrow',[nx2 nx1],[ny2 ny1],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
text(posM2_x_shifted + 0.02*diff(xlimv), posPeakVal/2, 'M_2', 'FontWeight','bold');

% --- M3: horizontal between pos and neg peaks (shifted slightly vertically) ---
[nx1, ny1] = to_norm(M3_x1, M3_y + shift_M3_y);
[nx2, ny2] = to_norm(M3_x2, M3_y + shift_M3_y);
annotation('line',[nx1 nx2],[ny1 ny2],'LineWidth',1.4);
annotation('arrow',[nx1 nx2],[ny1 ny2],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
annotation('arrow',[nx2 nx1],[ny2 ny1],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
text(mean([M3_x1 M3_x2]), M3_y + 2*shift_M3_y, 'M_3', 'HorizontalAlignment','center','FontWeight','bold');

% --- M4: vertical from 0 to negative peak, shifted right ---
if ~isempty(negPeakIdx)
    negM4_x_shifted = min(max(negPeakTime + shift_M4_x, xlimv(1)), xlimv(2));
    [nx1, ny1] = to_norm(negM4_x_shifted, 0);
    [nx2, ny2] = to_norm(negM4_x_shifted, negPeakVal);
    annotation('line',[nx1 nx2],[ny1 ny2],'LineWidth',1.4);
    annotation('arrow',[nx1 nx2],[ny1 ny2],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
    annotation('arrow',[nx2 nx1],[ny2 ny1],'HeadLength',6,'HeadWidth',6,'LineStyle','none');
    text(negM4_x_shifted + 0.02*diff(xlimv), negPeakVal/2, 'M_4', 'FontWeight','bold');
end

% final axis limits & aesthetics
xlim([min(t) max(t)]);
ylim(ylimv);  % keep original y-limits in case they were changed
drawnow;


%% Cluster-based permutation test (between groups, time domain)
% Inputs:
% - young_VEOG_all: [Nyoung x Ntime]
% - old_VEOG_all:   [Nold   x Ntime]
% - time:           [1 x Ntime]

ft_defaults;
rng(42);            % reproducibility

young_VEOG_all = young_out.VEOG_all;
old_VEOG_all   = old_out.VEOG_all;

Nyoung = size(young_VEOG_all,1);
Nold   = size(old_VEOG_all,1);

Ntotal = Nyoung + Nold;

time = time(:)';   % ensure row

% -------------------------------------------------------------------------
% Build FieldTrip timelock structures (single loop)
% -------------------------------------------------------------------------
all_VEOG = [young_VEOG_all; old_VEOG_all];   % [Ntotal x Ntime]
data = cell(1, Ntotal);

for s = 1:Ntotal
    data{s}.label  = {'VEOG'};
    data{s}.time   = time;
    data{s}.avg    = all_VEOG(s,:);           % [1 x time]
    data{s}.dimord = 'chan_time';
end

% -------------------------------------------------------------------------
% Design (independent samples)
% -------------------------------------------------------------------------
design = [ones(1,Nyoung), 2*ones(1,Nold)];  % 1 = young, 2 = old

% -------------------------------------------------------------------------
% Cluster-based permutation configuration
% -------------------------------------------------------------------------
cfg = [];
cfg.channel           = 'VEOG';
cfg.latency           = 'all';
cfg.method            = 'montecarlo';
cfg.statistic         = 'indepsamplesT';

cfg.correctm          = 'cluster';
cfg.clusteralpha      = 0.05;
cfg.clusterstatistic  = 'maxsum';

cfg.minnbchan         = 0;        % time-only clustering
cfg.neighbours        = [];

cfg.tail              = 0;        % two-sided
cfg.alpha             = 0.05;

cfg.numrandomization  = 10000;    % note: compute-heavy, 2k-5k often sufficient

cfg.design            = design;
cfg.ivar              = 1;        % group labels

% -------------------------------------------------------------------------
% Run statistics
% -------------------------------------------------------------------------
stat = ft_timelockstatistics(cfg, data{:});

% ============================
% Plot: Grand-average VEOG + SEM + cluster visualization
% ============================
figure('Color','w'); hold on;

% --- compute mean + SEM (use nan-aware std in case of NaNs) ---
mean_y = mean(young_VEOG_all,1);
sem_y  = nanstd(young_VEOG_all,0,1) ./ sqrt(max(1,size(young_VEOG_all,1)));

mean_o = mean(old_VEOG_all,1);
sem_o  = nanstd(old_VEOG_all,0,1) ./ sqrt(max(1,size(old_VEOG_all,1)));

% colors (match your lines below)
col_y = [1 0.2 0.2];   % red-ish
col_o = [0 0.2 0.8];   % blue-ish

% --- SEM patch coordinates ---
xp = [time, fliplr(time)];                 % closed polygon x
yp_y = [mean_y + sem_y, fliplr(mean_y - sem_y)];
yp_o = [mean_o + sem_o, fliplr(mean_o - sem_o)];

% Plot SEM bands first (so they sit at the very bottom)
semPatchYoung = patch(xp, yp_y, col_y, 'EdgeColor','none', 'FaceAlpha', 0.20);
semPatchOld   = patch(xp, yp_o, col_o, 'EdgeColor','none', 'FaceAlpha', 0.20);

% hide SEM patches from legend
set(semPatchYoung,'HandleVisibility','off');
set(semPatchOld,'HandleVisibility','off');

% ============================================================
% Overlay significant clusters (time-only) - same logic as yours
% create cluster patches *after* SEM so they appear above SEM
% ============================================================
legendPatch = gobjects(1,1);  % placeholder

if isfield(stat,'mask') && any(stat.mask(:))
    sigmask = stat.mask(1,:);
    d = diff([0 sigmask 0]);
    starts = find(d == 1);
    ends   = find(d == -1) - 1;
    yl = ylim;  % will get expanded automatically if needed

    sigTimeClusterPatch = gobjects(numel(starts),1);
    for k = 1:numel(starts)
        clusterT = mean(stat.stat(1, starts(k):ends(k)));
        if clusterT > 0
            clusterColor = [1 0.85 0.85];   % light red (young > old)
        else
            clusterColor = [0.85 0.85 1];   % light blue (old > young)
        end

        sigTimeClusterPatch(k) = patch( ...
            [time(starts(k)) time(ends(k)) time(ends(k)) time(starts(k))], ...
            [yl(1) yl(1) yl(2) yl(2)], ...
            clusterColor, ...
            'EdgeColor','none', ...
            'FaceAlpha',0.40);

        if k == 1
            legendPatch = sigTimeClusterPatch(k);
            legendPatch.DisplayName = 'Significant cluster';
        else
            sigTimeClusterPatch(k).HandleVisibility = 'off';
        end
    end
    % keep cluster patches visible (they are above SEM but below lines we will plot next)
end

% --- Now plot the group mean lines on top of the patches ---
youngLine = plot(time, mean_y, 'Color', col_y, 'LineWidth', 2);
oldLine   = plot(time, mean_o, 'Color', col_o, 'LineWidth', 2);

xlabel('Time (s)');
ylabel('Amplitude (µV)');
title('Grand-average VEOG');

% tidy up stacking: we want SEM (bottom) -> clusters (middle) -> lines (top)
axChildren = get(gca,'Children');  % current stacking order (top->bottom)
% Construct desired order: lines on top, then clusters, then sem patches, other items below
% Gather handles we created (some may be empty if clusters didn't exist)
linesH = [youngLine; oldLine];
clusterH = [];
if exist('sigTimeClusterPatch','var')
    clusterH = sigTimeClusterPatch(:);
end
semH = [semPatchYoung; semPatchOld];

% Build new child order filtering empty/invalid handles and keeping the rest of children after them
newOrder = [linesH(:); clusterH(:); semH(:)];
newOrder = newOrder(isgraphics(newOrder));   % remove non-graphics entries
% Append remaining children that we didn't explicitly order (so nothing disappears)
remaining = setdiff(axChildren, newOrder, 'stable');
set(gca,'Children',[newOrder; remaining]);

grid off;

% ---- Build legend once from handles ----
handles = [youngLine, oldLine];
labels  = {'Younger adults', 'Older adults'};

if ~isempty(legendPatch) && isgraphics(legendPatch)
    handles(end+1) = legendPatch;
    labels{end+1} = 'Significant cluster';
end

legend(handles, labels, 'Location','northwest', 'Box','off');