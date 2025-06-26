restoredefaultpath

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

% List of old subject IDs
old_subjects = { '8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45',	'5HO51',	'7MT51',	'4CH63',	'7CU61', '7WE49',	'6RE54',	'4KS63',	'6JH59',	'7HI40',	'6CR44',	'6CM48',	'10SH66',	'11MS53'};  % Replace with actual subject IDs
data_path = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging/Old';
out_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results';
% old {'8GH64', '8UH50', '7PE61',	'5MA56'	,'5CH43',	'7PS49', '6JU60', '5BN45',	'5HO51', '9RW37',	'7MT51',	'4CH63',	'7CU61', '7WE49',	'6RE54',	'4KS63',	'6JH59',	'7HI40',	'6CR44',	'6CM48',	'10SH66',	'11MS53'}
% exluded old '8GH64' '5CH43' '9RW37'

% Initialize storage for averaged results
old_freq_withERP = [];
old_freq_noERP = [];
old_mse_bin = [];

for iSub = 1:length(old_subjects)
    SUBJ = old_subjects{iSub};
    fprintf('Processing subject: %s\n', SUBJ);
    
    % Load cleaned data
    load(fullfile(data_path, [SUBJ '_cleaned_data.mat']), 'data_clean');
    
    % Rename the dataset to avoid confusion
    data_clean_old = data_clean;

    % Preserve trial information
    trl = data_clean_old.cfg.trl;

    % Step 1: Introduce TP7 as implicit reference
    cfg = [];
    cfg.channel = 'all';
    cfg.implicitref = 'TP7';
    cfg.reref = 'no';
    data_clean_old = ft_preprocessing(cfg, data_clean_old);
    
    % Step 2: Perform average referencing, ensuring TP7 remains
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = setdiff(data_clean_old.label, {'VEOG1', 'VEOG2', 'HEOG1', 'HEOG2'});
    cfg.channel = data_clean_old.label;
    data_clean_old = ft_preprocessing(cfg, data_clean_old);

    % Step 3: Reject trials with large variance
    disp 'reject trials with large variance' % last step, if avg ref messes trials up
    par = [];   par.badtrs = []; par.method = 'zscorecut';
    to_plot=0; if ispc; to_plot=1; end
    keeptrls = EM_ft_varcut3(data_clean_old, par, to_plot);
    cfg=[];
    cfg.trials = keeptrls;
    data_clean_old = ft_selectdata(cfg, data_clean_old);

    % Step 4: CSD transform (scalp current density)
    cfg = [];
    cfg.method = 'spline';
    cfg.elec = 'standard_1020.elc';
    data_clean_old = ft_scalpcurrentdensity(cfg, data_clean_old);

    % Step 5: Perform ERP analysis
    cfg = [];
    timelock = ft_timelockanalysis(cfg, data_clean_old);

    % Save ERP results
    erp_out_path = fullfile(out_path, 'timelock', SUBJ);
    mkdir(erp_out_path);
    save(fullfile(erp_out_path, 'timelock.mat'), 'timelock');

    % Create copy for ERP removal
    data_clean_old_noERP = data_clean_old;
    evoked = 'subtract';
    if strcmp(evoked, 'subtract')
        for itrial = 1:length(data_clean_old_noERP.trial)
            data_clean_old_noERP.trial{itrial} = data_clean_old_noERP.trial{itrial} - timelock.avg;
        end
    end

    % Step 6: Time-frequency analysis - With ERP
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'EEG';
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.foi = logspace(0, 1.5, 30); 
    cfg.t_ftimwin = ones(length(cfg.foi), 1) .* 0.5;
    cfg.toi = -1.25:0.05:1.25; 
    freq_withERP = ft_freqanalysis(cfg, data_clean_old);

    % Step 7: Time-frequency analysis - Without ERP
    freq_noERP = ft_freqanalysis(cfg, data_clean_old_noERP);

    % Step 8: Downsample data for entropy analysis
    cfg = [];
    cfg.resamplefs = 50;
    data_clean_old_noERP = ft_resampledata(cfg, data_clean_old_noERP);

    % Adjust trial information for downsampling
    downsample_factor = 350 / cfg.resamplefs;
    trl(:, 1:2) = round(trl(:, 1:2) / downsample_factor);
    data_clean_old_noERP.cfg.trl = trl;

    % Step 9: Entropy analysis
    addpath('C:\Users\morit\Desktop\Toolboxes_MATLAB\mMSE-master');
    cfg = [];
    cfg.m = 2;
    cfg.r = 0.5;
    cfg.timwin = 0.5;
    cfg.toi = -1.25:0.05:1.25;
    cfg.timescales = 1:8;
    cfg.recompute_r = 'perscale_toi_sp';
    cfg.coarsegrainmethod = 'filtskip';
    cfg.filtmethod = 'lp';
    cfg.mem_available = 20e+09;
    cfg.allowgpu = true;

    mse_bin = ft_entropyanalysis(cfg, data_clean_old_noERP);
    mse_bin.dimord = 'chan_freq_time'
    
    % Save individual results
    save(fullfile(out_path, [SUBJ '_freq_withERP.mat']), 'freq_withERP');
    save(fullfile(out_path, [SUBJ '_freq_noERP.mat']), 'freq_noERP');
    save(fullfile(out_path, [SUBJ '_mse_bin.mat']), 'mse_bin');

    % Accumulate results for averaging
    if isempty(old_freq_withERP)
        old_freq_withERP = freq_withERP;
        old_freq_noERP = freq_noERP;
        old_mse_bin = mse_bin;
    else
        old_freq_withERP.powspctrm = old_freq_withERP.powspctrm + freq_withERP.powspctrm;
        old_freq_noERP.powspctrm = old_freq_noERP.powspctrm + freq_noERP.powspctrm;
        old_mse_bin.sampen = old_mse_bin.sampen + mse_bin.sampen;
    end
end

% Average results across subjects
num_subjects = length(old_subjects);
old_freq_withERP.powspctrm = old_freq_withERP.powspctrm / num_subjects;
old_freq_noERP.powspctrm = old_freq_noERP.powspctrm / num_subjects;
old_mse_bin.sampen = old_mse_bin.sampen / num_subjects;

% Save averaged results
save(fullfile(out_path, 'old_average_freq_withERP.mat'), 'old_freq_withERP');
save(fullfile(out_path, 'old_average_freq_noERP.mat'), 'old_freq_noERP');
save(fullfile(out_path, 'old_average_mse_bin.mat'), 'old_mse_bin');

% Visualization of averaged results
cfg = [];
cfg.zlim = 'maxabs';
cfg.layout = 'EEG1010';
cfg.baseline = [-1 -0.5];
cfg.baselinetype = 'relchange';
cfg.showlabels = 'yes';
cfg.colorbar = 'yes';

% Time-frequency with ERP
figure;
ft_multiplotTFR(cfg, old_freq_withERP);
title('Average Time-Frequency Analysis with ERP');

% Time-frequency without ERP
figure;
ft_multiplotTFR(cfg, old_freq_noERP);
title('Average Time-Frequency Analysis without ERP');

% Entropy analysis
cfg.parameter = 'sampen';
cfg.ylim = [20 160];
cfg.xlim = [-1.5 1.5];


old_mse_bin.freq = old_mse_bin.timescales; % Map timescales to freq for compatibility

cfg = [];
cfg.parameter = 'sampen'; % Use the correct parameter for entropy visualization
cfg.layout = 'EEG1010';
cfg.zlim = 'maxabs';
cfg.colorbar = 'yes';
cfg.xlim = [-1.5 1.5];
cfg.ylim = [min(old_mse_bin.freq), max(old_mse_bin.freq)];

figure;
ft_multiplotTFR(cfg, old_mse_bin);
title('Average Entropy Analysis');


% List of young subject IDs
young_subjects = { '5IE98', '7KI87', '4ML96', '5SR93','11VZ96','10SH95', '6AU94','7FN98','8RS89','6LA93','10PM95', '8ST94','7SU95','7SN94', '6FM97', '9JN97','7JO97', '10PM95','11JI91','6KF96'};  % Replace with actual subject IDs
data_path = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging/Young';
out_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results';
% included young '5IE98', '7KI87', '4ML96', '5SR93','11VZ96','10SH95', '6AU94','7FN98','8RS89','6LA93','10PM95', '8ST94','7SU95','7SN94', '6FM97', '9JN97','7JO97', '10PM95','11JI91','6KF96'
% excluded young '5PH93' '10CH91', '6MS89'

% Initialize storage for averaged results
young_freq_withERP = [];
young_freq_noERP = [];
young_mse_bin = [];

for iSub = 1:length(young_subjects)
    SUBJ = young_subjects{iSub};
    fprintf('Processing subject: %s\n', SUBJ);
    
    % Load cleaned data
    load(fullfile(data_path, [SUBJ '_cleaned_data.mat']), 'data_clean');
    
    % Rename the dataset to avoid confusion
    data_clean_young = data_clean;

    % Preserve trial information
    trl = data_clean_young.cfg.trl;

    % Step 1: Introduce TP7 as implicit reference
    cfg = [];
    cfg.channel = 'all';
    cfg.implicitref = 'TP7';
    cfg.reref = 'no';
    data_clean_young = ft_preprocessing(cfg, data_clean_young);

    % Step 2: Perform average referencing, ensuring TP7 remains
    cfg = [];
    cfg.reref = 'yes';
    cfg.refchannel = setdiff(data_clean_young.label, {'VEOG1', 'VEOG2', 'HEOG1', 'HEOG2'});
    cfg.channel = data_clean_young.label;
    data_clean_young = ft_preprocessing(cfg, data_clean_young);

    % Step 3: Reject trials with large variance
    disp 'reject trials with large variance' % last step, if avg ref messes trials up
    par = [];   par.badtrs = []; par.method = 'zscorecut';
    to_plot=0; if ispc; to_plot=1; end
    keeptrls = EM_ft_varcut3(data_clean_young, par, to_plot);
    cfg=[];
    cfg.trials = keeptrls;
    data_clean_young = ft_selectdata(cfg, data_clean_young);

    % Step 4: CSD transform (scalp current density)
    cfg = [];
    cfg.method = 'spline';
    cfg.elec = 'standard_1020.elc';
    data_clean_young = ft_scalpcurrentdensity(cfg, data_clean_young);

    % Step 5: Perform ERP analysis
    cfg = [];
    timelock = ft_timelockanalysis(cfg, data_clean_young);

    % Save ERP results
    erp_out_path = fullfile(out_path, 'timelock', SUBJ);
    mkdir(erp_out_path);
    save(fullfile(erp_out_path, 'timelock.mat'), 'timelock');

    % Step 6: Create copy for ERP removal
    data_clean_young_noERP = data_clean_young;
    evoked = 'subtract';
    if strcmp(evoked, 'subtract')
        for itrial = 1:length(data_clean_young_noERP.trial)
            data_clean_young_noERP.trial{itrial} = data_clean_young_noERP.trial{itrial} - timelock.avg;
        end
    end

    % Step 7: Time-frequency analysis - With ERP
    cfg = [];
    cfg.output = 'pow';
    cfg.channel = 'EEG';
    cfg.method = 'mtmconvol';
    cfg.taper = 'hanning';
    cfg.foi = logspace(0, 1.5, 30); 
    cfg.t_ftimwin = ones(length(cfg.foi), 1) .* 0.5;
    cfg.toi = -1.25:0.05:1.25; 
    freq_withERP = ft_freqanalysis(cfg, data_clean_young);

    % Step 8: Time-frequency analysis - Without ERP
    freq_noERP = ft_freqanalysis(cfg, data_clean_young_noERP);

    % Step 9: Downsample data for entropy analysis
    cfg = [];
    cfg.resamplefs = 50;
    data_clean_young_noERP = ft_resampledata(cfg, data_clean_young_noERP);

    % Adjust trial information for downsampling
    downsample_factor = 350 / cfg.resamplefs;
    trl(:, 1:2) = round(trl(:, 1:2) / downsample_factor);
    data_clean_young_noERP.cfg.trl = trl;

    % Step 10: Entropy analysis
    addpath('C:\Users\morit\Desktop\Toolboxes_MATLAB\mMSE-master');
    cfg = [];
    cfg.m = 2;          % pattern parameter
    cfg.r = 0.5;        % similarity parameter
    cfg.timwin = 0.5;
    cfg.toi = -1.25:0.05:1.25;
    cfg.timescales = 1:8;
    cfg.recompute_r = 'perscale_toi_sp';
    cfg.coarsegrainmethod = 'filtskip';
    cfg.filtmethod = 'lp';
    cfg.mem_available = 20e+09;
    cfg.allowgpu = true;

    mse_bin = ft_entropyanalysis(cfg, data_clean_young_noERP);
    mse_bin.dimord = 'chan_freq_time'
    
    % Save individual results
    save(fullfile(out_path, [SUBJ '_freq_withERP.mat']), 'freq_withERP');
    save(fullfile(out_path, [SUBJ '_freq_noERP.mat']), 'freq_noERP');
    save(fullfile(out_path, [SUBJ '_mse_bin.mat']), 'mse_bin');

    % Accumulate results for averaging
    if isempty(young_freq_withERP)
        young_freq_withERP = freq_withERP;
        young_freq_noERP = freq_noERP;
        young_mse_bin = mse_bin;
    else
        young_freq_withERP.powspctrm = young_freq_withERP.powspctrm + freq_withERP.powspctrm;
        young_freq_noERP.powspctrm = young_freq_noERP.powspctrm + freq_noERP.powspctrm;
        young_mse_bin.sampen = young_mse_bin.sampen + mse_bin.sampen;
    end
end

% Average results across subjects
num_subjects = length(young_subjects);
young_freq_withERP.powspctrm = young_freq_withERP.powspctrm / num_subjects;
young_freq_noERP.powspctrm = young_freq_noERP.powspctrm / num_subjects;
young_mse_bin.sampen = young_mse_bin.sampen / num_subjects;

% Save averaged results
save(fullfile(out_path, 'young_average_freq_withERP.mat'), 'young_freq_withERP');
save(fullfile(out_path, 'young_average_freq_noERP.mat'), 'young_freq_noERP');
save(fullfile(out_path, 'young_average_mse_bin.mat'), 'young_mse_bin');

% Visualization of averaged results
cfg = [];
cfg.zlim = 'maxabs';
cfg.layout = 'EEG1010';
cfg.baseline = [-1 -0.5];
cfg.baselinetype = 'relchange';
cfg.showlabels = 'yes';
cfg.colorbar = 'yes';

% Time-frequency with ERP
figure;
ft_multiplotTFR(cfg, young_freq_withERP);
title('Average Time-Frequency Analysis with ERP');

% Time-frequency without ERP
figure;
ft_multiplotTFR(cfg, young_freq_noERP);
title('Average Time-Frequency Analysis without ERP');

% Entropy analysis
cfg.parameter = 'sampen';
cfg.ylim = [20 160];
cfg.xlim = [-1.5 1.5];


young_mse_bin.freq = young_mse_bin.timescales; % Map timescales to freq for compatibility

cfg = [];
cfg.parameter = 'sampen'; % Use the correct parameter for entropy visualization
cfg.layout = 'EEG1010';
cfg.zlim = 'maxabs';
cfg.colorbar = 'yes';
cfg.xlim = [-1.5 1.5];
cfg.ylim = [min(young_mse_bin.freq), max(young_mse_bin.freq)];

figure;
ft_multiplotTFR(cfg, young_mse_bin);
title('Average Entropy Analysis');


% Substract OA from YA TFA

addpath('C:\Users\morit\Desktop\FoPra_Daten\Analysis_Results');
addpath('C:\Users\morit\Desktop\Toolboxes_MATLAB\mMSE-master');
data_path = 'C:\Users\morit\Desktop\FoPra_Daten\Analysis_Results';

% Load the data files
load(fullfile(data_path, 'young_average_freq_noERP.mat'));
load(fullfile(data_path, 'old_average_freq_noERP.mat'));

% Step 2: Verify Structure Compatibility
% Ensure the data has the required fields for subtraction and visualization
disp(young_freq_noERP);
disp(old_freq_noERP);

% Check if dimord matches for subtraction
assert(strcmp(young_freq_noERP.dimord, old_freq_noERP.dimord), 'Dimension order mismatch!');
assert(isequal(size(young_freq_noERP.powspctrm), size(old_freq_noERP.powspctrm)), 'Data sizes mismatch!');

% Step 3: Compute Difference (Young - Old)
freq_diff = young_freq_noERP;
freq_diff.powspctrm = young_freq_noERP.powspctrm - old_freq_noERP.powspctrm; % Subtraction

% Step 4: Prepare Data for Visualization
cfg = [];
cfg.parameter = 'powspctrm';  % Power spectrum field
cfg.layout = 'EEG1010';       % Use the EEG1010 layout for multiplot
cfg.zlim = 'maxabs';          % Symmetrical color scaling
%cfg.baseline = [-1 -0.5];          Useful for Alpha suppression
%cfg.baselinetype = 'relchange';    Useful for Alpha Supression
cfg.colorbar = 'yes';         % Display colorbar
cfg.xlim = [-1.5 1.5];        % Time range
cfg.ylim = [min(young_freq_noERP.freq), max(young_freq_noERP.freq)]; % Frequency range

% Step 5: Visualize Young Group Time-Frequency Analysis
figure;
ft_multiplotTFR(cfg, young_freq_noERP);
title('Young Group Time-Frequency Analysis (No ERP)');

% Step 6: Visualize Old Group Time-Frequency Analysis
figure;
ft_multiplotTFR(cfg, old_freq_noERP);
title('Old Group Time-Frequency Analysis (No ERP)');

% Step 7: Visualize Time-Frequency Difference (Young - Old)
figure;
ft_multiplotTFR(cfg, freq_diff);
title('Time-Frequency Analysis Difference (Young - Old)');

% Substract OA from YA Entropy

% Step 1: Add Path and Load Data
addpath('C:\Users\morit\Desktop\FoPra_Daten\Analysis_Results');
data_path = 'C:\Users\morit\Desktop\FoPra_Daten\Analysis_Results';

% Load young and old average MSE data
load(fullfile(data_path, 'young_average_mse_bin.mat')); % Loads 'young_mse_bin'
load(fullfile(data_path, 'old_average_mse_bin.mat'));   % Loads 'old_mse_bin'

% Step 2: Prepare Data for Visualization
% Map timescales to freq for compatibility with FieldTrip plotting functions
young_mse_bin.freq = young_mse_bin.timescales;
old_mse_bin.freq = old_mse_bin.timescales;

% Step 3: Visualize Young Group Entropy
cfg = [];
cfg.parameter = 'sampen';       % Use the correct field for entropy visualization
cfg.layout = 'EEG1010';         % Use EEG1010 layout for multiplot
cfg.zlim = 'maxabs';            % Symmetrical color scaling
cfg.colorbar = 'yes';           % Display colorbar
cfg.xlim = [-1.5 1.5];          % Time range
cfg.ylim = [min(young_mse_bin.freq), max(young_mse_bin.freq)]; % Timescale range

figure;
ft_multiplotTFR(cfg, young_mse_bin);
title('Young Group Average Entropy');

% Step 4: Visualize Old Group Entropy
figure;
ft_multiplotTFR(cfg, old_mse_bin);
title('Old Group Average Entropy');

% Step 5: Compute and Visualize Entropy Difference (Young - Old)
mse_diff = young_mse_bin;
mse_diff.sampen = young_mse_bin.sampen - old_mse_bin.sampen; % Compute difference
mse_diff.freq = young_mse_bin.freq;                          % Ensure timescales are set
mse_diff.dimord = 'chan_freq_time';                          % Ensure correct dimension order

figure;
ft_multiplotTFR(cfg, mse_diff);
title('Entropy Difference (Young - Old)');