function eA_preproc_noblink(cfg)

%% preprocessing function
% assign the value of cfg.dataset to the variable dataset
dataset = cfg.dataset;
PREOUT_FIG = cfg.PREOUT_FIG;
PREOUT_CLEAN = cfg.PREOUT_CLEAN;
SUBJ = cfg.SUBJ;

cfg=[];                                     % initialize empty configuration structure
cfg.dataset = dataset;                      % previously defined dataset variable assigned to the dataset field of the cfg structure
cfg.channel = 'all';                        % sets channel field of the cfg structure to all, use all channels
cfg.continuous = 'yes';                     % specifies that data is continuous

data = ft_preprocessing(cfg);               % calls FieldTrip function using cfg as argument, perform preprocessing on data and save as such using the cfg structure 
                                            % look at the dataset using ft_data_browser(cfg,data)

%% Step 1: Linenoise removal Toolbox by Kluge

disp 'zapline-plus-main'                    % display string "" in the command window
cfg=[];                                     % initialize empty configuration structure
cfg.resample = 'yes';                       % data should be resampled
cfg.resamplefs = 350;                       % data sampling rate changed to 350Hz
cfg.detrend = 'no';                         % no detrending (removes trends and slow variations in the data)
data = ft_resampledata(cfg, data);          % uses the function resample on the data and saves this as data (change in data due to analysis)
dat = data.trial{1}(:,:);                   % select all rows from beginning to the third-to-last row; dat is now a numeric matrix: channels x time points

% call the zapline-function
[zaplineData, zaplineConfig, analyticsResults, plothandles] = clean_data_with_zapline_plus(dat, data.fsample, 'noisefreqs', 50);    % custom function passing cleaning parameters through it
                                                                                                                                    % dat - EEG data matrix from earlier
                                                                                                                                    % data.fsample - original sampling rate (350Hz)
                                                                                                                                    % noisefreqs, 50 - specify 50Hz as noise frequency to remove
saveas(gcf, fullfile(PREOUT_FIG, 'Controlanalysis',  sprintf('%s_zapline-plus_50Hz.png', SUBJ)) )       % Saves figures in path of PREOUT in folder figures

data.trial{1} = zaplineData;                                                                    % replaces original noisy channels with cleaned data

ft_postamble previous   data                                                                    % this copies the data.cfg structure into the cfg.previous field. You can also use it for multiple inputs, or for "varargin"
ft_postamble history    data                                                                    % this adds the local cfg structure to the output data structure, i.e. dataout.cfg = cfg

plotit=1;                                                                   % plot it 1 = TRUE
if plotit                                                                   
    disp 'check if zapline did well'                                        % message in command window
    cfgfreq              = [];                                              % initializes empty cfg structure
    cfgfreq.output       = 'pow';                                           % output should be power
    %   cfgfreq.channel      = 'all'; still applied to all channels
    cfgfreq.method       = 'mtmfft';                                        % multi-taper method with Fourier transform
    cfgfreq.taper        = 'hanning';                                       % reduce spectral leakage and smooths data; hanning window taper
    cfgfreq.keeptrials   = 'no';                                            % Use the spectrum and not individual trials
    cfgfreq.foilim       = [0 min(256/2, 200)];                             % frequencies of interest 0Hz to half of 256Hz(128Hz) or 200Hz
    cfgfreq.pad='nextpow2';                                                 % Zero-padding to increase resolution
    tempfreq = ft_freqanalysis(cfgfreq, data);                              % perform frequency analysis using ft_function on data
    figure; semilogy(tempfreq.freq, mean(tempfreq.powspctrm(1:end-3,:)))    % plots the average power spectrum across all but last 3 channels on a semilog y-axis
    clear tempfreq                                                          % remove the variable tempfreq from the workspace
end


%% Step 2: High-pass filter
disp('high pass filter the data')                                           % message in command window
cfg=[];                                                                     % initialize cfg structure 
cfg.hpfilter = 'yes';                                                       % high pass filter should be applied
cfg.hpfreq =  0.5;                                                          % everything above 0.5Hz is passed, below 0.5Hz is filtered
cfg.hpfiltord = 4;                                                          % order of the filter, higher orders mean steeper cutoff, determines the slope
cfg.continuous = 'yes';                                                     % data should be treated as continuous
data = ft_preprocessing(cfg, data);                                         % apply high-pass filter on data using ft_preprocessing


%% Step 3: Rereference Vert EOG to each other
cfg = [];                                                                   % initialize emtpy cfg structure
cfg.channel    = {'VEOG1', 'VEOG2'};                                        % specify channels of interest
cfg.reref      = 'yes';                                                     % enables re-referencing
cfg.refchannel = 'VEOG1';                                                   % calculate the signal relative to a specific reference channel (in this case VEOG1)
data_eogvert      = ft_preprocessing(cfg, data);                            % apply rereferencing via ft_preprocessing and data_eogvert is the output
data_eogvert.label{2} = 'EOGV';                                             % name this new reference channel EOGV


%% Extract the EOGV
cfg = [];                                                                                                                      
cfg.channel = 'EOGV';
data_eogvert   = ft_preprocessing(cfg, data_eogvert);                       % extract the channel EOGV from the data and name it data_eogvert


%% Step 4: Epoch the data into 3s segments
% Assuming 350 Hz sampling rate
fs = data.fsample;
epoch_length = 3; % in seconds
samples_per_epoch = fs * epoch_length;

nSamples = size(data.time{1}, 2);
nChunks = floor(nSamples / samples_per_epoch);

trl = zeros(nChunks, 3);
for i = 1:nChunks
    begsample = (i-1)*samples_per_epoch + 1;
    endsample = i*samples_per_epoch;
    offset = 0;
    trl(i,:) = [begsample endsample offset];
end

cfg = [];
cfg.trl = trl;
data_epoch = ft_redefinetrial(cfg, data);  % now 3s trials


%% Step 5: Detect blink-events on EOGV
cfg = [];
cfg.continuous = 'yes';
cfg.artfctdef.eog.channel = 'EOGV';
cfg.artfctdef.eog.bpfilter = 'yes';
cfg.artfctdef.eog.bpfreq = [1 15];
cfg.artfctdef.eog.hilbert = 'yes';
cfg.artfctdef.eog.cutoff = 2;
cfg.artfctdef.eog.feedback = 'no';

[cfg, artifact_blink] = ft_artifact_eog(cfg, data_eogvert);


%% Reject trials that overlap with blink events

cfg = [];
cfg.artfctdef.reject = 'complete';  % remove whole trials
cfg.artfctdef.eog.artifact = artifact_blink;

data_clean_3s = ft_rejectartifact(cfg, data_epoch);


%% Visual check

cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data_clean_3s);


%% Muscle artifact removal on blink-cleaned 1s trials
disp('Detecting muscle artifacts...');

% Prepare dummy trl from current 1s trials for artifact tools
trl = zeros(length(data_clean_3s.trial), 3);
for i = 1:length(data_clean_3s.trial)
    trl(i,1) = data_clean_3s.sampleinfo(i,1);
    trl(i,2) = data_clean_3s.sampleinfo(i,2);
    trl(i,3) = 0;
end

cfg = [];
cfg.continuous = 'yes';  % Required by artifact detection
cfg.trl = trl;

% Define channels to use (exclude EOGs)
eeg_channels = setdiff(data_clean_3s.label, {'VEOG1', 'VEOG2', 'HEOG1', 'HEOG2'});
cfg.artfctdef.muscle.channel = eeg_channels;
cfg.artfctdef.muscle.inspect = eeg_channels;

cfg.artfctdef.muscle.bpfilter = 'yes';
cfg.artfctdef.muscle.bpfreq = [80 100];
cfg.artfctdef.muscle.bpfiltord = 4;
cfg.artfctdef.muscle.hilbert = 'yes';
cfg.artfctdef.muscle.threshold = 4;       % Z-score threshold
cfg.artfctdef.muscle.cutoff = 20;         % Absolute cutoff for detection
cfg.artfctdef.muscle.feedback = 'yes';
cfg.artfctdef.muscle.boxcar = 0.2;

% Detect muscle artifacts
[cfg, artifact_muscle] = ft_artifact_muscle(cfg, data_clean_3s);

% Reject trials with muscle activity
cfg = [];
cfg.artfctdef.reject = 'complete';
cfg.artfctdef.muscle.artifact = artifact_muscle;

data_clean_3s_nomuscle = ft_rejectartifact(cfg, data_clean_3s);  % <-- final cleaned data

% (Optional) Save or plot
cfg = [];
cfg.viewmode = 'vertical';
ft_databrowser(cfg, data_clean_3s_nomuscle);


%% Perform ICA on blink- and muscle-cleaned data
disp('Run ICA...');

cfg = [];
cfg.channel = setdiff(data_clean_3s_nomuscle.label, {'VEOG2', 'HEOG2'});
data_clean_3s_nomuscle = ft_selectdata(cfg, data_clean_3s_nomuscle);

cfg = [];
cfg.channel = 'all';                          % Include all channels
cfg.method = 'runica';                        % ICA algorithm
cfg.runica.stop = 0.00000014;                 % Convergence threshold
comp = ft_componentanalysis(cfg, data_clean_3s_nomuscle);  % ICA decomposition

%% Step 1: Visualize time course of components
cfg = [];
cfg.viewmode = 'vertical';
cfg.channel = 1:min(60, size(comp.label,1));  % Show first 60 components or fewer if less exist
ft_databrowser(cfg, comp);

%% Step 2: Plot topographies of components
cfg = [];
cfg.component = 1:min(60, size(comp.label,1));
cfg.layout = 'EEG1010';                       % Update if your layout differs
cfg.viewmode = 'component';
cfg.comment = 'no';
cfg.marker = 'off';
figure('units', 'normalized', 'outerposition', [0 0 1 1]);
ft_topoplotIC(cfg, comp);

%% Step 3: Manual component rejection prompt
manual_reject_input = input('Enter components to reject (e.g., [1 2 4]): ');
if ~isempty(manual_reject_input)
    manual_reject = manual_reject_input;
else
    manual_reject = [];  % No components to reject
end

%% Step 4: Reject specified components and reconstruct data
cfg = [];
cfg.component = manual_reject;
data_clean = ft_rejectcomponent(cfg, comp, data_clean_3s_nomuscle);

%% Step 5: Preserve trial info and save cleaned data
if isfield(data_clean_3s_nomuscle.cfg, 'trl')
    data_clean.cfg.trl = data_clean_3s_nomuscle.cfg.trl;
else
    warning('No trial information (trl) found in source data.');
end

output_filename = fullfile(PREOUT_CLEAN, sprintf('%s_cleaned_data.mat', SUBJ));
save(output_filename, 'data_clean', '-v7.3');
disp(['Cleaned data saved to: ' output_filename]);

