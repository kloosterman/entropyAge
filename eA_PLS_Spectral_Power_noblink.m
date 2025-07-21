% PLS
restoredefaultpath;

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

%% Step 1: Add path to the behavioral measures file
addpath('C:/Users/morit/Desktop/FoPra_Daten');
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/plscmd-main');

behav = readtable('FoPra_Behavioral_Measures.xlsx');

% Apply zscore to columns 7-12
behav{:, 7:12} = zscore(behav{:, 7:12});

% Apply zscore to columns 14-19
behav{:, 14:19} = zscore(behav{:, 14:19});

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

%% Step 2: Load individual Spectral Power data for YA
young_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Young';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89', '4ML96', '5SR93', '11VZ96', '10SH95'};

young_freq = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_data_freq_young.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_freq{i} = temp.data_freq_young;            % Assuming the variable is called mse_bin
end

%% Step 1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'EEG1010';
neighbours = ft_prepare_neighbours(cfg);

%% Step 2: Configure statistics
cfg = [];
cfg.latency = [-1.5 1.5];                   % Time window for analysis
cfg.statistic = 'ft_statfun_indepsamplesT'; % Independent samples T-test
cfg.correctm = 'cluster';                   % Cluster correction method
cfg.numrandomization = 1000;                % Number of randomizations for the Monte Carlo test
cfg.method = 'montecarlo';                  % Monte Carlo method for statistics
cfg.clusteralpha = 0.05;                    % Cluster alpha level
cfg.clusterstatistic = 'maxsum';            % Statistic to use for clusters
cfg.tail = 0;                               % Two-tailed test
cfg.clustertail = 0;                        % Two-tailed cluster correction
cfg.alpha = 0.025;                          % Significance threshold
cfg.spmversion = 'spm12';                   % SPM12 version
cfg.neighbours = neighbours;


%% Step 3: Convert individual MSE data to freq-compatible format
% for i = 1:length(young_freq)
%     mse = young_freq{i};
%     freq_data_young = [];
%     freq_data_young.label = mse.label;               % Channel labels
%     freq_data_young.time = mse.time;                 % Time points
%     freq_data_young.freq = mse.timescales;           % Frequency (entropy timescales)
%     freq_data_young.powspctrm = mse.sampen;          % Entropy data as power spectrum
%     freq_data_young.dimord = 'chan_freq_time';       % Specify the dimension order
%     young_freq{i} = freq_data_young;                 % Replace with freq-compatible structure
% end

young_freq_all = ft_appendfreq([],young_freq{:})

%% Step 4: Configure statistics
cfg = [];
cfg.layout = 'EEG1010'
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 1000;                         % Number of permutation
cfg.num_boot = 1000;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 8;                            % 0 is Pearson corr, 8 is Spearman
cfg.design = behav_young;
cfg.num_cond = 1;                           % Number of conditions
cfg.num_subj_lst = size(young_freq_all.powspctrm,1); % Number of subjects per condition

%% Step 5: Compute statistics
stat_freq_young = ft_freqstatistics(cfg, young_freq_all);

%% Topoplot 100Hz doesn´t work
% Find index for 100 Hz
f_idx = find(stat_freq_young.freq == 100);

% Create dummy TFR-like structure for topoplot
dummy = [];
dummy.label     = stat_freq_young.label;
dummy.freq      = stat_freq_young.freq(f_idx);
dummy.time      = 0;
dummy.dimord    = 'chan_freq_time';
dummy.powspctrm = stat_freq_young.stat(:, f_idx);  % should be [nChan x 1]
dummy.powspctrm = reshape(dummy.powspctrm, [length(dummy.label), 1, 1]);  % enforce shape [chan x freq x time]

cfg = [];
cfg.layout    = 'EEG1010';
cfg.colorbar  = 'yes';
cfg.zlim      = 'maxabs';
cfg.colormap  = jet;

figure;
ft_topoplotTFR(cfg, dummy);
title('Topomap of PLS weights at 100 Hz');


%% Channel x Frequency

figure;
imagesc(stat_freq_young.freq, 1:length(stat_freq_young.label), stat_freq_young.stat);
xlabel('Frequency (Hz)');
ylabel('Channels');
yticks(1:length(stat_freq_young.label));
yticklabels(stat_freq_young.label);
title('PLS Stat Values: Channels vs Frequencies');
set(gca, 'YDir', 'normal');
colormap(jet); colorbar;

%% Brainscore x Behavscore
figure;
scatter(stat_freq_young.behavscores, -1*stat_freq_young.brainscores, 60, 'filled', ...
    'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0]);
lsline; % regression line
xlabel('Behavioral Scores');
ylabel('Brain Scores');
r = corr(stat_freq_young.behavscores, stat_freq_young.brainscores, 'type', 'Spearman');
title(sprintf('Brain vs Behavior Scores (Spearman = %.2f)', r));
grid on; box on;




%% Step 2: Load individual Spectral Power data for old
old_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Old';

% Define abbreviations
old_ids = {'8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45', '5HO51', ...
           '7MT51', '4CH63', '7CU61', '7WE49', '6RE54', '4KS63', '6JH59', ...
           '7HI40', '6CR44', '6CM48', '10SH66', '11MS53'};

old_freq = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_data_freq_old.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    old_freq{i} = temp.data_freq_old;            % Assuming the variable is called mse_bin
end

% Convert individual MSE data to freq-compatible format
% for i = 1:length(old_freq)
%     mse = old_freq{i};
%     freq_data_old = [];
%     freq_data_old.label = mse.label;               % Channel labels
%     freq_data_old.time = mse.time;                 % Time points
%     freq_data_old.freq = mse.timescales;           % Frequency (entropy timescales)
%     freq_data_old.powspctrm = mse.sampen;          % Entropy data as power spectrum
%     freq_data_old.dimord = 'chan_freq_time';       % Specify the dimension order
%     old_freq{i} = freq_data_old;                  % Replace with freq-compatible structure
% end

old_freq_all = ft_appendfreq([],old_freq{:})


% Step 2: Configure statistics
cfg = [];
cfg.layout = 'EEG1010'
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 1000;                         % Number of permutation
cfg.num_boot = 1000;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 8;                            % 0 is Pearson corr, 8 is Spearman
cfg.design = behav_old;
cfg.num_cond = 1;                           % Number of conditions
cfg.num_subj_lst = size(old_freq_all.powspctrm,1); % Number of subjects per condition

% Step 3: Compute statistics
stat_freq_old = ft_freqstatistics(cfg, old_freq_all);



%% Channel x Frequency

figure;
imagesc(stat_freq_old.freq, 1:length(stat_freq_old.label), stat_freq_old.stat);
xlabel('Frequency (Hz)');
ylabel('Channels');
yticks(1:length(stat_freq_old.label));
yticklabels(stat_freq_old.label);
title('PLS Stat Values: Channels vs Frequencies');
set(gca, 'YDir', 'normal');
colormap(jet); colorbar;

%% Brainscore x Behavscore
figure;
scatter(-1*stat_freq_old.behavscores, -1*stat_freq_old.brainscores, 60, 'filled', ...
    'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0]);                        % multiplied with -1 for interpretation
lsline; % regression line
xlabel('Behavioral Scores');
ylabel('Brain Scores');
r = corr(stat_freq_old.behavscores, stat_freq_old.brainscores, 'type', 'Spearman');
title(sprintf('Brain vs Behavior Scores (Spearman = %.2f)', r));
grid on; box on;