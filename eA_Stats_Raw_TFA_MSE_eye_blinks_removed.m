%% Statistical analysis of power spectrum and MSE between YA and OA

% Reset MATLAB path to avoid conflicts
restoredefaultpath;

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

% Add Statistical Parametric Mapping Toolbox (SPM12)
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/spm12');

% Add Path and Load Data
addpath('C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results');
data_path = 'C:/Users/morit/Desktop/FoPra_Daten/Conntrolanalysis_Results';
young_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Young';
old_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Old';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '4ML96', '5SR93', '11VZ96', '10SH95'};

old_ids = {'8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45', '5HO51', ...
           '7MT51', '4CH63', '7CU61', '7WE49', '6RE54', '4KS63', '6JH59', ...
           '7HI40', '6CR44', '6CM48', '10SH66', '11MS53'};

% Load individual Frequency data for young and old
young_freq = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_data_freq_young.mat']);
    temp = load(file_path);                                                 % Adjust this if the variable in the file is named differently
    young_freq{i} = temp.data_freq_young;                                        % Assuming the variable is called cleaned_data
end

old_freq = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_data_freq_old.mat']);
    temp = load(file_path);                                                 % Adjust this if the variable in the file is named differently
    old_freq{i} = temp.data_freq_old;                                        % Assuming the variable is called cleaned_data
end

% Concatenate young_freq data along the 4th dimension (rpt)
cfg = [];
young_freq = ft_appendfreq(cfg, young_freq{:}); % Concatenate all subjects
young_freq.powspctrm = permute(young_freq.powspctrm, [1, 2, 3, 4]); % Ensure subjects are in 4th dimension

% Concatenate old_freq data along the 4th dimension (rpt)
cfg = [];
old_freq = ft_appendfreq(cfg, old_freq{:}); % Concatenate all subjects
old_freq.powspctrm = permute(old_freq.powspctrm, [1, 2, 3, 4]); % Ensure subjects are in 4th dimension

% Check the output
disp(young_freq);
disp(old_freq);

%% - - - Start the analysis - - - %%

%% Time-Frequency-Analysis - Power Analysis

% Step 1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'EEG1010';
neighbours = ft_prepare_neighbours(cfg);

% Step 2: Configure statistics
cfg = [];
cfg.latency = [0 1];                        % Time window for analysis
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

% Design matrix for the independent samples t-test
cfg.design = [ones(1, length(young_ids)), 2 * ones(1, length(old_ids))];    % Young = 1, Old = 2
cfg.ivar = 1;                                                               % The grouping factor is in the first column of the design

% Step 3: Compute statistics
stat_freq = ft_freqstatistics(cfg, young_freq, old_freq);

% Display cluster information
disp(stat_freq.posclusters);
disp(stat_freq.negclusters);

% Step 4: Configure cluster plot
cfg = [];
cfg.alpha = 0.025;                                          % Set the significance level for clusters
cfg.parameter = 'stat';                                     % Specify the parameter to plot (e.g., t-statistics)
cfg.zlim = 'maxabs';                                        % Set color limits for visualization
cfg.layout = 'EEG1010';                                     % Specify your layout (use the same as earlier)
cfg.marker = 'on';                                          % Show electrode markers
cfg.highlight = 'on';                                       % Highlight the significant channels
cfg.highlightchannel = stat_freq.posclusterslabelmat == 1;  % Only highlight channels in the first positive cluster
cfg.colorbar = 'yes';                                       % Include colorbar for scale reference

% Plot the positive cluster = higher power in younger adults
ft_multiplotTFR(cfg, stat_freq);

% Plot MultiplotTFR for Significant Clusters
cfg = [];
cfg.layout = 'EEG1010';                                                     % Use the layout specified earlier
cfg.parameter = 'stat';                                                     % Use the statistic (e.g., t-values)
cfg.colorbar = 'yes';                                                       % Include a colorbar
cfg.zlim = 'maxabs';                                                        % Adjust the color axis limits
stat_freq.mask = double(stat_freq.posclusterslabelmat == 1);                % Highlight positive cluster
cfg.maskparameter = 'mask';                                                 % Use the mask for plotting significant areas

ft_multiplotTFR(cfg, stat_freq);


%% Entropy analysis

% Load individual Entropy data for young and old

young_mse = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                                                 % Adjust this if the variable in the file is named differently
    young_mse{i} = temp.mse_bin;                                            % Assuming the variable is called mse_bin
end

old_mse = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                                                 % Adjust this if the variable in the file is named differently
    old_mse{i} = temp.mse_bin;                                              % Assuming the variable is called mse_bin
end

% Convert individual MSE data to freq-compatible format
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

% Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'EEG1010';                         % Use the layout you already confirmed
neighbours = ft_prepare_neighbours(cfg);

% Configure statistics
cfg = [];
cfg.latency = [0 1];                            % Adjusted to match your MSE time range
cfg.frequency = [20 110];                       % Frequency range for analysis (entropy timescales)
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

% Design matrix for the independent samples t-test
cfg.design = [ones(1, length(young_ids)), 2 * ones(1, length(old_ids))];    % Young = 1, Old = 2
cfg.ivar = 1;                                                               % The grouping factor is in the first column of the design

% Compute statistics
stat_mse = ft_freqstatistics(cfg, young_mse, old_mse);

% Display cluster information
disp(stat_mse.posclusters);
disp(stat_mse.negclusters);

% Configure cluster plot
cfg = [];
cfg.alpha = 0.025;                                                          % Set the significance level for clusters
cfg.parameter = 'stat';                                                     % Specify the parameter to plot (e.g., t-statistics)
cfg.zlim = 'maxabs';                                                        % Set color limits for visualization
cfg.layout = 'EEG1010';                                                     % Specify your layout (use the same as earlier)
cfg.marker = 'on';                                                          % Show electrode markers
cfg.highlight = 'on';                                                       % Highlight the significant channels
cfg.highlightchannel = find(stat_mse.negclusterslabelmat == 1);             % Channels in the first positive cluster
cfg.colorbar = 'yes';                                                       % Include colorbar for scale reference

ft_multiplotTFR(cfg,stat_mse)

% Plot MultiplotTFR for Significant Clusters
cfg = [];
cfg.layout = 'EEG1010';                                                     % Use the layout specified earlier
cfg.parameter = 'stat';                                                     % Use the statistic (e.g., t-values)
cfg.colorbar = 'yes';                                                       % Include a colorbar
cfg.zlim = 'maxabs';                                                        % Adjust the color axis limits
stat_mse.mask = double(stat_mse.negclusterslabelmat == 1);                  % Highlight negative cluster
cfg.maskparameter = 'mask';                                                 % Use the mask for plotting significant areas

ft_multiplotTFR(cfg, stat_mse);                                             % higher Entropy for OA all electrodes

% Clusterplot for Raw Power Analysis
figure
cfg=[];
cfg.layout = 'EEG1010';
cfg.clus2plot = 1;
load colormap_jetlightgray.mat
cfg.colormap = cmap
cfg.integratetype = 'trapz';                                    % mean or trapz
cfg.clussign = 'pos';
cfg.parameter = 'stat';
cfg.titleTFR = 'Power YA vs OA';
stat_mse.mask = double(stat_freq.posclusterslabelmat == 1);     % Highlight negative cluster
cfg.maskparameter = 'mask';                                     % Use the mask for plotting significant areas
ft_clusterplot3D(cfg, stat_freq);                               % Greater Power for YA compared to OA

figure
cfg=[];
cfg.layout = 'EEG1010';
cfg.clus2plot = 1;
load colormap_jetlightgray.mat
cfg.colormap = cmap
cfg.integratetype = 'trapz';                                    % mean or trapz
cfg.clussign = 'neg';
cfg.parameter = 'stat';
cfg.titleTFR = 'Power YA vs OA';
stat_mse.mask = double(stat_freq.negclusterslabelmat == 1);     % Highlight negative cluster
cfg.maskparameter = 'mask';                                     % Use the mask for plotting significant areas
ft_clusterplot3D(cfg, stat_freq);                               % Greater Power for OA compared to YA

% Clusterplot for MSE analysis
figure
cfg=[];
cfg.layout = 'EEG1010';
cfg.clus2plot = 1;
load colormap_jetlightgray.mat
cfg.colormap = cmap
cfg.integratetype = 'trapz';                                    % mean or trapz
cfg.clussign = 'neg';
cfg.parameter = 'stat';
cfg.titleTFR = 'mse YA vs OA';
stat_mse.mask = double(stat_mse.negclusterslabelmat == 1);      % Highlight negative cluster
cfg.maskparameter = 'mask';                                     % Use the mask for plotting significant areas
ft_clusterplot3D(cfg, stat_mse);                                % Greater MSE for OA compared to YA

figure
cfg=[];
cfg.layout = 'EEG1010';
cfg.clus2plot = 1;
load colormap_jetlightgray.mat
cfg.colormap = cmap
cfg.integratetype = 'trapz';                                    % mean or trapz
cfg.clussign = 'pos';
cfg.parameter = 'stat';
cfg.titleTFR = 'mse YA vs OA';
stat_mse.mask = double(stat_mse.posclusterslabelmat == 1);      % Highlight positive cluster
cfg.maskparameter = 'mask';                                     % Use the mask for plotting significant areas
ft_clusterplot3D(cfg, stat_mse);                                % greater MSE for YA compared to OA

%% Visualize YA Entropy and OA Entropy separately for every channel (scalp)
% Assuming `young_mse` and `old_mse` are already loaded in the workspace

% Create the configuration structure for plotting
cfg = [];
cfg.layout = 'EEG1010';             % Use the EEG1010 layout to display all channels (ensure it's available in your FieldTrip setup)
load colormap_jetlightgray.mat
cfg.colormap = cmap
cfg.channel = 'all';                % Plot all channels
cfg.frequency = [20 100];           % Frequency range (use the full range or a specific one)
cfg.time = young_mse.time;          % Use the time vector from your data
cfg.zlim = [0.8 1.1];               % Set color limits to the maximum absolute value (optional)
cfg.xlim = [0 1];                   % Time window
cfg.ylim = [20 100];                % Frequency range (from young_mse.freq)

% Prepare the frequency structure for young_mse
young_mse_freq = young_mse;
% If you have multiple repetitions, average over them (optional)
young_mse_freq.powspctrm = mean(young_mse.powspctrm, 1);    % Average across repetitions

% Plot using ft_multiplotTFR for young_mse
figure;
ft_multiplotTFR(cfg, young_mse_freq);                       % MSE across whole time of YA

% Do the same for old_mse
old_mse_freq = old_mse;
old_mse_freq.powspctrm = mean(old_mse.powspctrm, 1);        % Average across repetitions (optional)

% Plot using ft_multiplotTFR for old_mse
figure;
ft_multiplotTFR(cfg, old_mse_freq);                         % MSE across whole time of OA


%% What is that for?

% Step 1: Average across repetitions (first dimension)
young_mse_avg = young_mse;
young_mse_avg.powspctrm = mean(young_mse.powspctrm, 1); % Average across repetitions

old_mse_avg = old_mse;
old_mse_avg.powspctrm = mean(old_mse.powspctrm, 1); % Average across repetitions

% Step 2: Compute the difference (YA - OA)
diff_mse = young_mse_avg;
diff_mse.powspctrm = young_mse_avg.powspctrm - old_mse_avg.powspctrm;

% Step 3: Extract only Fz
channel_idx = find(strcmp(diff_mse.label, 'Fz')); % Find index of Fz
diff_mse_Fz = diff_mse;
diff_mse_Fz.powspctrm = diff_mse.powspctrm(:, channel_idx, :, :); % Keep only Fz data
diff_mse_Fz.label = {'Fz'};  % Update label

% Step 4: Plot using ft_singleplotTFR
cfg = [];
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.channel = 'Fz';          % Show only Fz
cfg.frequency = [20 100];     % Frequency range
cfg.xlim = [-1.5 1.5];       % Time window
cfg.ylim = [20 100];         % Frequency range (from young_mse.freq)
cfg.zlim = 'maxabs';         % Keep color scaling consistent


%% YA and OA Entropy plots Again channel Fz

% Create the configuration structure for plotting
cfg = [];
cfg.channel = 'Fz';                % Plot only the Fz channel
load colormap_jetlightgray.mat
cfg.colormap = cmap
cfg.frequency = [20 100];           % Frequency range (use the full range or a specific one)
cfg.time = young_mse.time;          % Use the time vector from your data
cfg.zlim = [0.8 1.1];               % Set color limits to the maximum absolute value (optional)
cfg.xlim = [0 1];                   % Time window
cfg.ylim = [20 100];                % Frequency range (from young_mse.freq)

% Prepare the frequency structure for young_mse
young_mse_freq = young_mse;
% If you have multiple repetitions, average over them (optional)
young_mse_freq.powspctrm = mean(young_mse.powspctrm, 1);  % Average across repetitions

% Plot using ft_singleplotTFR for young_mse
figure;
ft_singleplotTFR(cfg, young_mse_freq);
xlabel('Time (s)', 'FontSize', 10);  % Set x-label
ylabel('Time scale', 'FontSize', 10); % Set y-label
title('Young Adults Fz', 'FontSize', 12); % Add title for young adults

% Do the same for old_mse
old_mse_freq = old_mse;
old_mse_freq.powspctrm = mean(old_mse.powspctrm, 1);  % Average across repetitions (optional)

% Plot using ft_singleplotTFR for old_mse
figure;
ft_singleplotTFR(cfg, old_mse_freq);
xlabel('Time (s)', 'FontSize', 10);  % Set x-label
ylabel('Time scale', 'FontSize', 10); % Set y-label
title('Old Adults Fz', 'FontSize', 12); % Add title for old adults

% Create the configuration structure for plotting
cfg = [];
cfg.channel = 'Fz';                % Plot only the Fz channel
load colormap_jetlightgray.mat
cfg.colormap = cmap
cfg.frequency = [20 100];           % Frequency range (use the full range or a specific one)
cfg.time = young_mse.time;          % Use the time vector from your data
cfg.zlim = [0.8 1.1];              % Set color limits to the maximum absolute value (optional)
cfg.xlim = [-1.5 1.5];             % Time window
cfg.ylim = [20 100];               % Frequency range (from young_mse.freq)

% Prepare the frequency structure for young_mse
young_mse_freq = young_mse;
% If you have multiple repetitions, average over them (optional)
young_mse_freq.powspctrm = mean(young_mse.powspctrm, 1);  % Average across repetitions

% Prepare the frequency structure for old_mse
old_mse_freq = old_mse;
old_mse_freq.powspctrm = mean(old_mse.powspctrm, 1);  % Average across repetitions (optional)

