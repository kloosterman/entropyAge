%% PLS analysis MSE for comparing the conditions blink and noblink

restoredefaultpath;

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

%% Add path to the behavioral measures file
addpath('C:/Users/morit/Desktop/FoPra_Daten');
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/plscmd-main');

behav = readtable('FoPra_Behavioral_Measures.xlsx');

behavNames = behav.Properties.VariableNames;
behav_young = behav(behav.Age_Group==1,:);
behav_young = table2array(behav_young(:,7:19));
behavNames = behavNames(7:19);
behav_young = behav_young(:,[1:6, 9:end]);
behavNames = behavNames([1:6, 9:end]);

behavNames = behav.Properties.VariableNames;
behav_old = behav(behav.Age_Group==2,:);
behav_old = table2array(behav_old(:,7:19));
behavNames = behavNames(7:19);
behav_old = behav_old(:,[1:6, 9:end]);
behavNames = behavNames([1:6, 9:end]);

%% PLS analysis MSE for YA 
% Load individual Entropy data for YA blink data
young_path_blink = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Young';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89' '4ML96', '5SR93', '11VZ96', '10SH95'};

young_mse_blink = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path_blink, [young_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_mse_blink{i} = temp.young_mse_bin;            % Assuming the variable is called mse_bin
end

%Load individual Entropy data for YA noblink data
young_path_noblink = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Young';

young_mse_noblink = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path_noblink, [young_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_mse_noblink{i} = temp.mse_bin;            % Assuming the variable is called mse_bin
end

%% Step 1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'EEG1010';
neighbours = ft_prepare_neighbours(cfg);

% Convert individual MSE blink data to freq-compatible format
for i = 1:length(young_mse_blink)
    mse = young_mse_blink{i};
    freq_data_young = [];
    freq_data_young.label = mse.label;               % Channel labels
    freq_data_young.time = mse.time;                 % Time points
    freq_data_young.freq = mse.timescales;           % Frequency (entropy timescales)
    freq_data_young.powspctrm = mse.sampen;          % Entropy data as power spectrum
    freq_data_young.dimord = 'chan_freq_time';       % Specify the dimension order
    young_mse_blink{i} = freq_data_young;                  % Replace with freq-compatible structure
end

young_mse_all_blink = ft_appendfreq([],young_mse_blink{:})

%% Use specific timing where the effect of stats is the highest
%% Extract single timepoint (-0.4s) from blink data
time_of_interest = -0.4;

% Find closest time index
[~, idx_time] = min(abs(young_mse_all_blink.time - time_of_interest));

% Extract powspctrm at that time
powspctrm_toi = young_mse_all_blink.powspctrm(:,:,:,idx_time);  % [trials x chan x freq]

% Create new struct matching noblink format
young_mse_all_blink_toi = struct();
young_mse_all_blink_toi.label = young_mse_all_blink.label;
young_mse_all_blink_toi.time  = young_mse_all_blink.time(idx_time);  % scalar
young_mse_all_blink_toi.freq  = young_mse_all_blink.freq;
young_mse_all_blink_toi.dimord = 'rpt_chan_freq_time';
young_mse_all_blink_toi.powspctrm = powspctrm_toi;
young_mse_all_blink_toi.cfg = young_mse_all_blink.cfg;

% Convert individual MSE noblink data to freq-compatible format
for i = 1:length(young_mse_noblink)
    mse = young_mse_noblink{i};
    freq_data_young = [];
    freq_data_young.label = mse.label;               % Channel labels
    freq_data_young.time = mse.time;                 % Time points
    freq_data_young.freq = mse.timescales;           % Frequency (entropy timescales)
    freq_data_young.powspctrm = mse.sampen;          % Entropy data as power spectrum
    freq_data_young.dimord = 'chan_freq_time';       % Specify the dimension order
    young_mse_noblink{i} = freq_data_young;                  % Replace with freq-compatible structure
end

young_mse_all_noblink = ft_appendfreq([],young_mse_noblink{:})

% Number of subjects per condition
nSubBlink   = size(young_mse_all_blink_toi.powspctrm,1);
nSubNoBlink = size(young_mse_all_noblink.powspctrm,1);

% % Flatten features: [channels × timescales] → 1D per subject
% features = size(young_mse_all_blink_toi.powspctrm,2) * size(young_mse_all_blink_toi.powspctrm,3);
% 
% datamat_lst = cell(1,2);
% % Each cell: [subjects × features] now
% datamat_lst{1} = reshape(young_mse_all_blink_toi.powspctrm, nSubBlink, features);  
% datamat_lst{2} = reshape(young_mse_all_noblink.powspctrm, nSubNoBlink, features);

% Design vector: 1 = blink, 2 = noblink
design_matrix = [ones(1, nSubBlink), 2*ones(1, nSubNoBlink)];

young_mse_all_concat = young_mse_all_blink_toi;  % start with blink data
young_mse_all_concat.powspctrm = cat(1, young_mse_all_blink_toi.powspctrm, young_mse_all_noblink.powspctrm);  % concatenate along subjects
young_mse_all_concat.time = young_mse_all_blink_toi.time;  % keep single timepoint

% Run PLS
cfg = [];
cfg.layout = 'EEG1010';
cfg.statistic = 'ft_statfun_pls';
cfg.num_perm = 1000;
cfg.num_boot = 1000;
cfg.method = 'analytic';
cfg.pls_method = 1;               % taskPLS
cfg.design = design_matrix
cfg.num_cond = 2;
cfg.num_subj_lst = [nSubBlink, nSubNoBlink];

stat_mse_young_condition = ft_statfun_pls(cfg, young_mse_all_concat);



disp('Blink data size:')
size(young_mse_all_blink_toi.powspctrm)

disp('No-blink data size:')
size(young_mse_all_noblink.powspctrm)

disp('Concatenated data size:')
size(young_mse_all_concat.powspctrm)

disp('Design matrix length:')
length(design_matrix)

disp('num_subj_lst:')
[nSubBlink, nSubNoBlink]