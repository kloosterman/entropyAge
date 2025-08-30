%% PLS Interaction analysis for Entropy in blink data
restoredefaultpath;

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

%% Step 1: Add path to the behavioral measures file
addpath('C:/Users/morit/Desktop/FoPra_Daten');
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/plscmd-main');

behav = readtable('FoPra_Behavioral_Measures.xlsx');

behavNames = behav.Properties.VariableNames;
behav_young = behav(behav.Age_Group==1,:);
behav_young = table2array(behav_young(:,7:19));
behavNames = behavNames(7:19);
behav_young = behav_young(:,[1:6, 9:end]);
behavNames = behavNames([1:6, 9:end]);

behav_young_demeaned = behav_young - mean(behav_young, 1);                  % demeaning is important for the interaction analysis

behavNames = behav.Properties.VariableNames;
behav_old = behav(behav.Age_Group==2,:);
behav_old = table2array(behav_old(:,7:19));
behavNames = behavNames(7:19);
behav_old = behav_old(:,[1:6, 9:end]);
behavNames = behavNames([1:6, 9:end]);

behav_old_demeaned = behav_old - mean(behav_old, 1);                        % demeaning is an important for the interaction analysis

behav_all_demeaned = [behav_young_demeaned; behav_old_demeaned];

%% Step 2: Load individual Entropy data for YA
young_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Young';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89', '4ML96', '5SR93', '11VZ96', '10SH95'};

young_mse = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_mse{i} = temp.young_mse_bin;            % Assuming the variable is called old_mse_bin
end

%% Step 3: Convert individual MSE data to freq-compatible format
for i = 1:length(young_mse)
    mse = young_mse{i};
    freq_data_young = [];
    freq_data_young.label = mse.label;               % Channel labels
    freq_data_young.time = mse.time;                 % Time points
    freq_data_young.freq = mse.timescales;           % Frequency (entropy timescales)
    freq_data_young.powspctrm = mse.sampen;          % Entropy data as power spectrum
    freq_data_young.dimord = 'chan_freq_time';       % Specify the dimension order
    young_mse{i} = freq_data_young;                  % Replace with freq-compatible structure
end

young_mse_all = ft_appendfreq([],young_mse{:})

% Create a new copy to preserve the original
young_mse_all_demeaned = young_mse_all;

% Extract and demean powspctrm across subjects (dimension 1)
mean_pow_young = mean(young_mse_all.powspctrm, 1);  % [1 x 60 x 8]
young_mse_all_demeaned.powspctrm = young_mse_all.powspctrm - mean_pow_young;


%% Step 2: Load individual Entropy data for old
old_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Old';

% Define abbreviations
old_ids = {'8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45', '5HO51', ...
           '7MT51', '4CH63', '7CU61', '7WE49', '6RE54', '4KS63', '6JH59', ...
           '7HI40', '6CR44', '6CM48', '10SH66', '11MS53'};

old_mse = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    old_mse{i} = temp.mse_bin;            % Assuming the variable is called mse_bin
end

%% Step 3: Convert individual MSE data to freq-compatible format
for i = 1:length(old_mse)
    mse = old_mse{i};
    freq_data_old = [];
    freq_data_old.label = mse.label;               % Channel labels
    freq_data_old.time = mse.time;                 % Time points
    freq_data_old.freq = mse.timescales;           % Frequency (entropy timescales)
    freq_data_old.powspctrm = mse.sampen;          % Entropy data as power spectrum
    freq_data_old.dimord = 'chan_freq_time';       % Specify the dimension order
    old_mse{i} = freq_data_old;                  % Replace with freq-compatible structure
end

old_mse_all = ft_appendfreq([],old_mse{:})

% Create a new copy to preserve the original
old_mse_all_demeaned = old_mse_all;

% Extract and demean powspctrm across subjects (dimension 1)
mean_pow_old = mean(old_mse_all.powspctrm, 1);  % [1 x 60 x 8]
old_mse_all_demeaned.powspctrm = old_mse_all.powspctrm - mean_pow_old;

%% Step 4: Calculate group PLS analysis (Interaction)

cfg = [];
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';                               % PLS statistics
cfg.num_perm = 1000;                                            % Number of permutation
cfg.num_boot = 1000;
cfg.method = 'analytic';                                        % analytic method for statistics
cfg.pls_method = 3;                                             % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 8;                                                % 0 is Pearson corr, 8 is Spearman
cfg.num_cond = 1;                                               % Number of conditions
cfg.interaction = 'yes';                                        % add group interaction to the model
cfg.contrast = [-1 1];                                          % -1 for YA, 1 for OA
cfg.design = behav_all_demeaned                                 % append behav OA to YA
cfg.num_subj_lst = [20 19];                                     % Number of subjects per condition, array!
stat_mse_2group = ft_freqstatistics(cfg, young_mse_all_demeaned, old_mse_all_demeaned);

%% Step 5: Plot Interaction results

%% Step 5.1: Plot Topoplot colapsed over timescale
cfg = [];
cfg.layout = 'EEG1010';
cfg.parameter = 'stat';
cfg.interactive = 'no';
cfg.zlim = 'maxabs';
stat_mse_2group.mask = double(stat_mse_2group.stat > 3 | stat_mse_2group.stat < -3);
load colormap_jetlightgray.mat
cfg.colormap = cmap;
%cfg.maskparameter = 'mask';
cfg.colorbar = 'yes';
ft_topoplotER(cfg, stat_mse_2group);

%% Step 5.2: Plot Entropy x Timescale
figure;
plot(stat_mse_2group.freq, stat_mse_2group.stat);
xlabel('Frequency (Hz)');
ylabel('Statistical Value');
title('Statistical Values over Frequency');
grid on;