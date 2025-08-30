%% Cluster-based Permutation Test (cbPT) and other analysis of MSE blink data
restoredefaultpath;

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

% Add Statistical Parametric Mapping Toolbox (SPM12)
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/spm12');

% Add Path and Load Data
addpath('C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results');
young_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Young';
old_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Old';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89', '4ML96', '5SR93', '11VZ96', '10SH95'};

old_ids = {'8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45', '5HO51', ...
           '7MT51', '4CH63', '7CU61', '7WE49', '6RE54', '4KS63', '6JH59', ...
           '7HI40', '6CR44', '6CM48', '10SH66', '11MS53'};

%% Step 1: Load individual Entropy data for young and old

young_mse = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                                                 % Adjust this if the variable in the file is named differently
    young_mse{i} = temp.young_mse_bin;                                            % Assuming the variable is called mse_bin
end

old_mse = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                                                 % Adjust this if the variable in the file is named differently
    old_mse{i} = temp.old_mse_bin;                                              % Assuming the variable is called mse_bin
end

%% Step 2: Convert individual MSE data to freq-compatible format
for i = 1:length(young_mse)
    mse = young_mse{i};
    freq_data = [];
    freq_data.label = mse.label;               % Channel labels
    freq_data.time = mse.time;                 % Time points
    freq_data.freq = mse.timescales;           % Frequency (entropy timescales)
    freq_data.powspctrm = mse.sampen;          % Entropy data as power spectrum
    freq_data.dimord = 'chan_freq_time';       % Specify the dimension order
    young_mse{i} = freq_data;                  % Replace with freq-compatible structure
end

for i = 1:length(old_mse)
    mse = old_mse{i};
    freq_data = [];
    freq_data.label = mse.label;               % Channel labels
    freq_data.time = mse.time;                 % Time points
    freq_data.freq = mse.timescales;           % Frequency (entropy timescales)
    freq_data.powspctrm = mse.sampen;          % Entropy data as power spectrum
    freq_data.dimord = 'chan_freq_time';       % Specify the dimension order
    old_mse{i} = freq_data;                    % Replace with freq-compatible structure
end

% Concatenate young_mse data along the 4th dimension (subjects)
cfg = [];
young_mse = ft_appendfreq(cfg, young_mse{:});

% Concatenate old_mse data along the 4th dimension (subjects)
cfg = [];
old_mse = ft_appendfreq(cfg, old_mse{:});

% Check the output
disp(young_mse);
disp(old_mse);

%% Step 3: Cluster based permutation test

% Step 3.1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'Custom_EEG1010.mat';                         % Use the layout you already confirmed
neighbours = ft_prepare_neighbours(cfg);

% Step 3.2: Configure statistics for cluster based permutation test
cfg = [];
cfg.latency = [-1.5 1.5];                       % Adjusted to match your MSE time range
cfg.frequency = [20 120];                       % Frequency range for analysis (entropy timescales)
cfg.statistic = 'ft_statfun_indepsamplesT';     % Independent samples T-test
cfg.correctm = 'cluster';                       % Cluster correction method
cfg.numrandomization = 1000;                    % Number of randomizations for the Monte Carlo test
cfg.method = 'montecarlo';                      % Monte Carlo method for statistics
cfg.clusteralpha = 0.05;                        % Cluster alpha level
cfg.clusterstatistic = 'maxsum';                % Statistic to use for clusters
cfg.tail = 0;                                   % Two-tailed test
cfg.clustertail = 0;                            % Two-tailed cluster correction
cfg.alpha = 0.025;                              % Significance threshold
cfg.spmversion = 'spm12';                       % SPM12 version
cfg.neighbours = neighbours;

% Step 3.3: Design matrix for the independent samples t-test
cfg.design = [ones(1, length(young_ids)), 2 * ones(1, length(old_ids))];    % Young = 1, Old = 2
cfg.ivar = 1;                                                               % The grouping factor is in the first column of the design

% Step 3.4: Compute statistics
stat_mse = ft_freqstatistics(cfg, young_mse, old_mse);

%% Step 4: Plotting the results

% Step 4.1: Plot positive (higher entropy for YA) and negative (higher entropy for OA) clusters for whole scalp
cfg = [];
cfg.layout = 'Custom_EEG1010.mat';                                                     % Use the layout specified earlier
cfg.alpha = 0.025;
cfg.parameter = 'stat';                                                     % Use the statistic (e.g., t-values)
cfg.colorbar = 'yes';                                                       % Include a colorbar
cfg.zlim = 'maxabs';                                                        % Adjust the color axis limits
cfg.maskparameter = 'mask';  
load colormap_jetlightgray.mat
cfg.colormap = cmap

stat_mse.mask = double(stat_mse.negclusterslabelmat == 1);                  % Highlight negative cluster                                                % Use the mask for plotting significant areas
ft_multiplotTFR(cfg, stat_mse);                                             % higher Entropy for OA all electrodes

stat_mse.mask = double(stat_mse.posclusterslabelmat == 1);
ft_multiplotTFR(cfg, stat_mse);


% Step 4.2: Plot positive and negative cluster over all electrodes
figure
cfg=[];
cfg.layout = 'Custom_EEG1010.mat';
cfg.clus2plot = 1;
cfg.colormap = cmap
cfg.integratetype = 'trapz';                                    % mean or trapz
cfg.clussign = 'neg';
cfg.parameter = 'stat';
cfg.titleTFR = 'MSE Negative Cluster';
cfg.maskparameter = 'mask';                                     % Use the mask for plotting significant areas

stat_mse.mask = double(stat_mse.negclusterslabelmat == 1);      % Highlight negative cluster
ft_clusterplot3D(cfg, stat_mse);                                % Greater MSE for OA compared to YA

figure
cfg=[];
cfg.layout = 'Custom_EEG1010.mat';
cfg.clus2plot = 1;
cfg.colormap = cmap
cfg.intergratetype = 'trapz';
cfg.clussign = 'pos';
cfg.parameter = 'stat';
cfg.titleTFR = 'MSE Positive Cluster'
stat_mse.mask = double(stat_mse.posclusterslabelmat == 1);
ft_clusterplot3D(cfg, stat_mse);