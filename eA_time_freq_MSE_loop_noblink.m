%% Time-frequency-analysis for YA and OA separately (looped)

% 
restoredefaultpath

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

%% List of old subject IDs
old_subjects = { '8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45',	'5HO51',	'7MT51',	'4CH63',	'7CU61', '7WE49',	'6RE54',	'4KS63',	'6JH59',	'7HI40',	'6CR44',	'6CM48',	'10SH66',	'11MS53'};  % Replace with actual subject IDs
data_path = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging_Controlanalysis/Old';
out_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Old';

% old {'8GH64', '8UH50', '7PE61',	'5MA56'	,'5CH43',	'7PS49', '6JU60', '5BN45',	'5HO51', '9RW37',	'7MT51',	'4CH63',	'7CU61', '7WE49',	'6RE54',	'4KS63',	'6JH59',	'7HI40',	'6CR44',	'6CM48',	'10SH66',	'11MS53'}
% exluded old '8GH64' '5CH43' '9RW37'

% Initialize storage for averaged results
avg_old_freq = [];
old_mse_bin = [];

for iSub = 1:length(old_subjects)
    SUBJ = old_subjects{iSub};
    fprintf('Processing subject: %s\n', SUBJ);
    
    % Load cleaned data
    load(fullfile(data_path, [SUBJ '_cleaned_data.mat']), 'data_clean');
    
    % Rename the dataset to avoid confusion
    data_clean_old = data_clean;

    %% Step 1: Introduce TP7 as implicit reference
    cfg = [];
    cfg.channel = 'all';
    cfg.implicitref = 'TP7';
    cfg.reref = 'no';
    data_clean_old = ft_preprocessing(cfg, data_clean_old);

    %% Step 2: Perform average referencing, ensuring TP7 remains
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = setdiff(data_clean_old.label, {'HEOG1'});
    cfg.channel = data_clean_old.label;
    data_clean_old = ft_preprocessing(cfg, data_clean_old);

    %% Step 3: Reject trials with large variance
    disp 'reject trials with large variance'                                % last step, if avg ref messes trials up
    par = [];   par.badtrs = []; par.method = 'zscorecut';
    to_plot=0; if ispc; to_plot=1; end
    keeptrls = EM_ft_varcut3(data_clean_old, par, to_plot);
    cfg=[];
    cfg.trials = keeptrls;
    data_clean_old = ft_selectdata(cfg, data_clean_old);

    % Preserve trial information
    trl = data_clean_old.sampleinfo;

    %% Step 4: CSD transform (scalp current density)
    cfg = [];
    cfg.method = 'spline';
    cfg.elec = 'standard_1020.elc';
    data_clean_old = ft_scalpcurrentdensity(cfg, data_clean_old);

    %% Step 5: Limit trial count for both TFR + MSE to make it comparable to the blink-event analysis
    max_trials = 200;
    nTrials = length(data_clean_old.trial);
    if nTrials > max_trials
    cfg = [];
    cfg.trials = randperm(nTrials, max_trials);
    data_clean_old = ft_selectdata(cfg, data_clean_old);
    fprintf('Trial count reduced from %d to %d for subject %s\n', nTrials, max_trials, SUBJ);
    end
 
    %% Step 6: Time-frequency analysis
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'EEG';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi        = 2:2:100;                           % 2 to 100Hz 
    cfg.t_ftimwin  = 5 ./ cfg.foi;                      % window length
    cfg.t_ftimwin(cfg.t_ftimwin > 1) = 1;               % Limit window size
    cfg.toi        = 0:0.05:1;                          % centers of windows that fully fit inside the trial -> when using 0:0.5:1 you have to add padding; 0.25:0.05:0.75 works without padding
    cfg.pad        = 2;                                 % ← add this line to pad trial to 2 seconds
    cfg.keeptrials = 'yes';

    data_freq_old = ft_freqanalysis(cfg, data_clean_old);

    %% Step 7: Downsample data for entropy analysis
    cfg = [];
    cfg.resamplefs = 50;
    data_MSE_old = ft_resampledata(cfg, data_clean_old);

    % Adjust trial information for downsampling
    downsample_factor = 350 / cfg.resamplefs;
    trl(:, 1:2) = round(trl(:, 1:2) / downsample_factor);
    data_MSE_old.cfg.trl = trl;

    % if isempty(data_MSE_old.trial)
    % warning('No trials after resampling for subject %s', SUBJ);
    % continue;
    % end

    %% Step 8: Entropy analysis
    addpath('C:\Users\morit\Desktop\Toolboxes_MATLAB\mMSE-master');
    cfg = [];
    cfg.m = 2;
    cfg.r = 0.5;
    cfg.timwin = 0.5;
    cfg.toi = 0.25:0.05:0.75;
    cfg.timescales = 1:8;
    cfg.recompute_r = 'perscale_toi_sp';
    cfg.coarsegrainmethod = 'filtskip';
    cfg.filtmethod = 'lp';
    cfg.mem_available = 20e+09;
    cfg.allowgpu = true;

    mse_bin = ft_entropyanalysis(cfg, data_MSE_old);
    mse_bin.dimord = 'chan_freq_time'

    % Convert to subject average
    cfg = [];
    cfg.keeptrials = 'no';
    data_freq_old = ft_freqdescriptives(cfg, data_freq_old);

    % Save individual results
    save(fullfile(out_path, [SUBJ '_data_freq_old.mat']), 'data_freq_old');
    save(fullfile(out_path, [SUBJ '_mse_bin.mat']), 'mse_bin');

end

% Initialize cell array
old_freq_all = cell(1, length(old_subjects));

% Load all freq_clean files
old_freq_all = {};
for i = 1:length(old_subjects)
    subj = old_subjects{i};
    file = fullfile(out_path, [subj '_data_freq_old.mat']);
    if exist(file, 'file')
        tmp = load(file);
        if isfield(tmp, 'data_freq_old') && isstruct(tmp.data_freq_old)
            old_freq_all{end+1} = tmp.data_freq_old;
        else
            warning('No valid data_freq_old in %s', file);
        end
    else
        warning('File missing: %s', file);
    end
end

% Average over subjects
cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';                                                  % set to 'yes' if you want to preserve subject dimension
avg_old_freq = ft_freqgrandaverage(cfg, old_freq_all{:});

% Average over time
cfg = [];
cfg.avgovertime = 'yes';
cfg.nanmean = 'yes';                                                        % <--- IMPORTANT: ignores NaNs in time-avg
avg_old_freq = ft_selectdata(cfg, avg_old_freq);

% Average over channels
avg_old_power = mean(avg_old_freq.powspctrm, 1, 'omitnan');                 % 1 x nFreq average over all channels


figure;
plot(avg_old_freq.freq, avg_old_power, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Average Power Spectrum (2–100 Hz)');
xlim([0 100]);        % <-- this ensures full frequency range is shown
grid on;


% Plot the results TFR
% cfg = [];
% cfg.layout       = 'EEG1010.lay';     % Adjust if you use a different system
% cfg.zlim         = 'maxabs';          % or [0 5] if you want fixed power scale
% cfg.baseline     = [0.25 0.35];       % Optional: define baseline (match your toi)
% cfg.baselinetype = 'absolute';        % or 'absolute', 'db', 'relchange'
% cfg.showlabels   = 'yes';             % Show channel labels
% 
% ft_multiplotTFR(cfg, avg_old_freq);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize cell array
old_mse_all = cell(1, length(old_subjects));

for i = 1:length(old_subjects)
    subj = old_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % loads variable 'mse_bin'
    
    mse = tmp.mse_bin;   % shortcut
    
    % Rename 'sampen' to 'powspctrm'
    mse.powspctrm = mse.sampen;
    
    % Assign 'freq' field (FieldTrip expects this)
    mse.freq = mse.timescales;
    
    % Optionally remove 'sampen' to avoid confusion
    mse = rmfield(mse, 'sampen');
    
    % Make sure dimord is correct (should already be)
    mse.dimord = 'chan_freq_time';
    
    % Save back to the cell array
    old_mse_all{i} = mse;
end


% Initialize cell array
old_mse_all = cell(1, length(old_subjects));

for i = 1:length(old_subjects)
    subj = old_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % should contain variable 'mse_bin'
    old_mse_all{i} = tmp.mse_bin;
end

% Convert to FieldTrip-compatible freq-like structure
for i = 1:length(old_mse_all)
    old_mse_all{i}.powspctrm = old_mse_all{i}.sampen;
    old_mse_all{i}.freq = old_mse_all{i}.timescales;    % Important!
    old_mse_all{i}.dimord = 'chan_freq_time';
    old_mse_all{i} = rmfield(old_mse_all{i}, 'sampen'); % Remove original
end

% Grand average across subjects
cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';
avg_old_mse = ft_freqgrandaverage(cfg, old_mse_all{:});

% Average entropy over channels and time (ignoring NaNs)
gavg_old_mse = squeeze(mean(mean(avg_old_mse.powspctrm, 1, 'omitnan'), 3, 'omitnan'));

% Timescales on x-axis
x = avg_old_mse.freq;

% Plot
figure;
plot(x, gavg_old_mse, 'LineWidth', 2);
xlabel('Timescale');
ylabel('Entropy (Sample Entropy)');
title('Average Entropy across Timescales');
grid on;




%% List of young subject IDs
young_subjects = { '5IE98', '6AU94', '7FN98', '8RS89','6LA93','10PM95', '8ST94','7SU95','7SN94', '6FM97', '11JI91', '9JN97','7JO97', '6KF96', '7KI87','6MS89' ,'4ML96', '5SR93','11VZ96','10SH95'};  % Replace with actual subject IDs
data_path = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging_Controlanalysis/Young';
out_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Young';
% included young '5IE98', '7KI87', '4ML96', '5SR93','11VZ96','10SH95', '6AU94','7FN98','8RS89','6LA93','10PM95', '8ST94','7SU95','7SN94', '6FM97', '9JN97','7JO97', '10PM95','11JI91','6KF96'
% excluded young '5PH93' '10CH91'

% Initialize storage for averaged results
avg_young_freq = [];
young_mse_bin = [];

for iSub = 1:length(young_subjects)
    SUBJ = young_subjects{iSub};
    fprintf('Processing subject: %s\n', SUBJ);
    
    % Load cleaned data
    load(fullfile(data_path, [SUBJ '_cleaned_data.mat']), 'data_clean');
    
    % Rename the dataset to avoid confusion
    data_clean_young = data_clean;

    %% Step 1: Introduce TP7 as implicit reference
    cfg = [];
    cfg.channel = 'all';
    cfg.implicitref = 'TP7';
    cfg.reref = 'no';
    data_clean_young = ft_preprocessing(cfg, data_clean_young);

    %% Step 2: Perform average referencing, ensuring TP7 remains
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = setdiff(data_clean_young.label, {'HEOG1'});
    cfg.channel = data_clean_young.label;
    data_clean_young = ft_preprocessing(cfg, data_clean_young);

    %% Step 3: Reject trials with large variance
    disp 'reject trials with large variance'                                % last step, if avg ref messes trials up
    par = [];   par.badtrs = []; par.method = 'zscorecut';
    to_plot=0; if ispc; to_plot=1; end
    keeptrls = EM_ft_varcut3(data_clean_young, par, to_plot);
    cfg=[];
    cfg.trials = keeptrls;
    data_clean_young = ft_selectdata(cfg, data_clean_young);

    % Preserve trial information
    trl = data_clean_young.sampleinfo;

    %% Step 4: CSD transform (scalp current density)
    cfg = [];
    cfg.method = 'spline';
    cfg.elec = 'standard_1020.elc';
    data_clean_young = ft_scalpcurrentdensity(cfg, data_clean_young);

    %% Step 5: Limit trial count for both TFR + MSE to make it comparable to the blink-event analysis
    max_trials = 100;
    nTrials = length(data_clean_young.trial);
    if nTrials > max_trials
    cfg = [];
    cfg.trials = randperm(nTrials, max_trials);
    data_clean_young = ft_selectdata(cfg, data_clean_young);
    fprintf('Trial count reduced from %d to %d for subject %s\n', nTrials, max_trials, SUBJ);
    end

    %% Step 6: Time-frequency analysis
    cfg = [];
    cfg.output     = 'pow';
    cfg.channel    = 'EEG';
    cfg.method     = 'mtmconvol';
    cfg.taper      = 'hanning';
    cfg.foi        = 2:2:100;                   % 2 to 100Hz 
    cfg.t_ftimwin  = 5 ./ cfg.foi;              % window length
    cfg.t_ftimwin(cfg.t_ftimwin > 1) = 1;       % Limit window size
    cfg.toi        = 0:0.05:1;                  % centers of windows that fully fit inside the trial -> when using 0:0.5:1 you have to add padding; 0.25:0.05:0.75 works without padding
    cfg.pad        = 2;                         % ← add this line to pad trial to 2 seconds
    cfg.keeptrials = 'yes';

    data_freq_young = ft_freqanalysis(cfg, data_clean_young);

    %% Step 7: Downsample data for entropy analysis
    cfg = [];
    cfg.resamplefs = 50;
    data_MSE_young = ft_resampledata(cfg, data_clean_young);

    % Adjust trial information for downsampling
    downsample_factor = 350 / cfg.resamplefs;
    trl(:, 1:2) = round(trl(:, 1:2) / downsample_factor);
    data_MSE_young.cfg.trl = trl;

    % if isempty(data_MSE_young.trial)
    % warning('No trials after resampling for subject %s', SUBJ);
    % continue;
    % end

    %% Step 8: Entropy analysis
    addpath('C:\Users\morit\Desktop\Toolboxes_MATLAB\mMSE-master');
    cfg = [];
    cfg.m = 2;
    cfg.r = 0.5;
    cfg.timwin = 0.5;
    cfg.toi = 0.25:0.05:0.75;
    cfg.timescales = 1:8;
    cfg.recompute_r = 'perscale_toi_sp';
    cfg.coarsegrainmethod = 'filtskip';
    cfg.filtmethod = 'lp';
    cfg.mem_available = 20e+09;
    cfg.allowgpu = true;

    mse_bin = ft_entropyanalysis(cfg, data_MSE_young);
    mse_bin.dimord = 'chan_freq_time'

    % Convert to subject average
    cfg = [];
    cfg.keeptrials = 'no';
    data_freq_young = ft_freqdescriptives(cfg, data_freq_young);

    % Save individual results
    save(fullfile(out_path, [SUBJ '_data_freq_young.mat']), 'data_freq_young');
    save(fullfile(out_path, [SUBJ '_mse_bin.mat']), 'mse_bin');

end

% Initialize cell array
young_freq_all = cell(1, length(young_subjects));

% Load all freq_clean files
young_freq_all = {};
for i = 1:length(young_subjects)
    subj = young_subjects{i};
    file = fullfile(out_path, [subj '_data_freq_young.mat']);
    if exist(file, 'file')
        tmp = load(file);
        if isfield(tmp, 'data_freq_young') && isstruct(tmp.data_freq_young)
            young_freq_all{end+1} = tmp.data_freq_young;
        else
            warning('No valid data_freq_young in %s', file);
        end
    else
        warning('File missing: %s', file);
    end
end

% Average over subjects
cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';                                                  % set to 'yes' if you want to preserve subject dimension
avg_young_freq = ft_freqgrandaverage(cfg, young_freq_all{:});

% Average over time
cfg = [];
cfg.avgovertime = 'yes';
cfg.nanmean = 'yes';                                                        % <--- IMPORTANT: ignores NaNs in time-avg
avg_young_freq = ft_selectdata(cfg, avg_young_freq);

% Average over channels
avg_young_power = mean(avg_young_freq.powspctrm, 1, 'omitnan');                 % 1 x nFreq average over all channels


figure;
plot(avg_young_freq.freq, avg_young_power, 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power');
title('Average Power Spectrum (2–100 Hz)');
xlim([0 100]);                                                              % <-- this ensures full frequency range is shown
grid on;


% Initialize cell array
young_mse_all = cell(1, length(young_subjects));

for i = 1:length(young_subjects)
    subj = young_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % loads variable 'mse_bin'
    
    mse = tmp.mse_bin;   % shortcut
    
    % Rename 'sampen' to 'powspctrm'
    mse.powspctrm = mse.sampen;
    
    % Assign 'freq' field (FieldTrip expects this)
    mse.freq = mse.timescales;
    
    % Optionally remove 'sampen' to avoid confusion
    mse = rmfield(mse, 'sampen');
    
    % Make sure dimord is correct (should already be)
    mse.dimord = 'chan_freq_time';
    
    % Save back to the cell array
    young_mse_all{i} = mse;
end


% Initialize cell array
young_mse_all = cell(1, length(young_subjects));

for i = 1:length(young_subjects)
    subj = young_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % should contain variable 'mse_bin'
    young_mse_all{i} = tmp.mse_bin;
end

% Convert to FieldTrip-compatible freq-like structure
for i = 1:length(young_mse_all)
    young_mse_all{i}.powspctrm = young_mse_all{i}.sampen;
    young_mse_all{i}.freq = young_mse_all{i}.timescales;    % Important!
    young_mse_all{i}.dimord = 'chan_freq_time';
    young_mse_all{i} = rmfield(young_mse_all{i}, 'sampen'); % Remove original
end

% Grand average across subjects
cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';
avg_young_mse = ft_freqgrandaverage(cfg, young_mse_all{:});

% Average entropy over channels and time (ignoring NaNs)
gavg_young_mse = squeeze(mean(mean(avg_young_mse.powspctrm, 1, 'omitnan'), 3, 'omitnan'));

% Timescales on x-axis
x = avg_young_mse.freq;

% Plot
figure;
plot(x, gavg_young_mse, 'LineWidth', 2);
xlabel('Timescale');
ylabel('Entropy (Sample Entropy)');
title('Average Entropy across Timescales');
grid on;