function out = extract_ant_post_timecourses(cfg, stat, dataYA, dataOA)
% EXTRACT_ANT_POST_TIMECOURSES
% Extract anterior and posterior time courses from chan x freq x time data.
%
% INPUT
%   cfg.layout          = FieldTrip layout name or layout struct
%   cfg.mask            = logical mask, size chan x freq x time
%                         If empty, all channels contribute equally.
%   cfg.timerange       = [tmin tmax] for defining channel weights
%                         default = full time range
%   cfg.freqrange       = [fmin fmax] for defining channel weights
%                         default = full freq range
%   cfg.weighting       = 'cluster' or 'equal' (default = 'cluster')
%   cfg.splitmode       = 'medianY' or 'labels' (default = 'medianY')
%   cfg.ant_labels      = cellstr, only used if splitmode='labels'
%   cfg.post_labels     = cellstr, only used if splitmode='labels'
%   cfg.normalizeweights = true/false (default = true)
%
%   stat                = struct with stat.label, stat.freq, stat.time
%   dataYA              = YA data, size subj x chan x freq x time
%   dataOA              = OA data, size subj x chan x freq x time
%
% OUTPUT
%   out.ant.YA          = YA anterior time courses, subj x time
%   out.post.YA         = YA posterior time courses, subj x time
%   out.ant.OA          = OA anterior time courses, subj x time
%   out.post.OA         = OA posterior time courses, subj x time
%   out.weights         = channel weights used
%   out.anterior_idx    = logical channel index
%   out.posterior_idx   = logical channel index
%   out.chanY           = y positions from layout
%
% NOTES
% - Time courses are averaged over selected freqrange and channels.
% - If weighting='cluster', weights come from cluster occupancy in the
%   selected freq/time window.
% - Same ROIs are used for YA and OA.

% ---------------- options ----------------
layout_in          = ft_getopt(cfg, 'layout', []);
mask               = ft_getopt(cfg, 'mask', []);
timerange          = ft_getopt(cfg, 'timerange', [min(stat.time) max(stat.time)]);
freqrange          = ft_getopt(cfg, 'freqrange', [min(stat.freq) max(stat.freq)]);
weighting          = ft_getopt(cfg, 'weighting', 'cluster');
splitmode          = ft_getopt(cfg, 'splitmode', 'medianY');
ant_labels         = ft_getopt(cfg, 'ant_labels', {});
post_labels        = ft_getopt(cfg, 'post_labels', {});
normalizeweights   = ft_getopt(cfg, 'normalizeweights', true);

% ---------------- checks ----------------
if isempty(layout_in)
    error('cfg.layout is required');
end
if ~isfield(stat, 'label') || ~isfield(stat, 'freq') || ~isfield(stat, 'time')
    error('stat must contain label, freq, and time');
end
if ndims(dataYA) ~= 4 || ndims(dataOA) ~= 4
    error('dataYA and dataOA must be subj x chan x freq x time');
end

nChan = numel(stat.label);

% ---------------- layout ----------------
tmpcfg = [];
tmpcfg.layout = layout_in;
layout = ft_prepare_layout(tmpcfg, stat);

[selchan, sellay] = match_str(stat.label, layout.label);
if numel(selchan) ~= nChan
    error('Could not match all stat.label entries to layout.label');
end

chanY = nan(nChan,1);
chanY(selchan) = layout.pos(sellay,2);

% ---------------- define anterior/posterior ----------------
switch lower(splitmode)
    case 'mediany'
        ymed = median(chanY, 'omitnan');
        anterior_idx  = chanY > ymed;
        posterior_idx = chanY < ymed;

    case 'labels'
        anterior_idx  = ismember(stat.label, ant_labels);
        posterior_idx = ismember(stat.label, post_labels);

    otherwise
        error('Unknown cfg.splitmode: %s', splitmode);
end

% remove channels without layout position from both sets
anterior_idx(isnan(chanY))  = false;
posterior_idx(isnan(chanY)) = false;

if ~any(anterior_idx)
    error('No anterior channels selected');
end
if ~any(posterior_idx)
    error('No posterior channels selected');
end

% ---------------- window selection ----------------
fsel = stat.freq >= freqrange(1) & stat.freq <= freqrange(2);
tsel = stat.time >= timerange(1) & stat.time <= timerange(2);

if ~any(fsel)
    error('No freq bins in requested freqrange');
end
if ~any(tsel)
    error('No time bins in requested timerange');
end

% ---------------- channel weights ----------------
switch lower(weighting)
    case 'cluster'
        if isempty(mask)
            warning('cfg.weighting=''cluster'' but no mask provided; using equal weights');
            weights = ones(nChan,1);
        else
            if ~isequal(size(mask), [nChan numel(stat.freq) numel(stat.time)])
                error('mask size must be chan x freq x time');
            end
            tmpmask = mask(:, fsel, tsel);
            weights = squeeze(sum(sum(double(tmpmask),3),2));
        end

    case 'equal'
        weights = ones(nChan,1);

    otherwise
        error('Unknown cfg.weighting: %s', weighting);
end

weights = weights(:);

% keep only weights within each region later
if normalizeweights
    if sum(weights(anterior_idx)) > 0
        weights_ant = weights;
        weights_ant(~anterior_idx) = 0;
        weights_ant = weights_ant ./ sum(weights_ant(anterior_idx));
    else
        weights_ant = double(anterior_idx);
        weights_ant = weights_ant ./ sum(weights_ant);
    end

    if sum(weights(posterior_idx)) > 0
        weights_post = weights;
        weights_post(~posterior_idx) = 0;
        weights_post = weights_post ./ sum(weights_post(posterior_idx));
    else
        weights_post = double(posterior_idx);
        weights_post = weights_post ./ sum(weights_post);
    end
else
    weights_ant = double(anterior_idx);
    weights_post = double(posterior_idx);
end

% ---------------- extract time courses ----------------
out.ant.YA  = local_extract_tc(dataYA, fsel, weights_ant);
out.post.YA = local_extract_tc(dataYA, fsel, weights_post);

out.ant.OA  = local_extract_tc(dataOA, fsel, weights_ant);
out.post.OA = local_extract_tc(dataOA, fsel, weights_post);

% ---------------- extras ----------------
out.weights        = weights;
out.weights_ant    = weights_ant;
out.weights_post   = weights_post;
out.anterior_idx   = anterior_idx;
out.posterior_idx  = posterior_idx;
out.chanY          = chanY;
out.time           = stat.time;
out.freqsel        = fsel;
out.timesel        = tsel;
out.label          = stat.label;

end

% ========================================================================
function tc = local_extract_tc(data, fsel, weights_chan)
% data: subj x chan x freq x time
% weights_chan: chan x 1

nSubj = size(data,1);
nTime = size(data,4);

% average over selected frequencies first
dat_f = squeeze(mean(data(:,:,fsel,:), 3));   % subj x chan x time

% weighted channel average for each time point
tc = nan(nSubj, nTime);

for isub = 1:nSubj
    x = squeeze(dat_f(isub,:,:));  % chan x time
    tc(isub,:) = weights_chan' * x;
end
end