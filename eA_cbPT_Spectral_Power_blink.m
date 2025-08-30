%% Cluster-based Permutation Test (cbPT) and other analysis of Spectral Power blink data

% Reset MATLAB path to avoid conflicts
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

% Load individual Frequency data for young and old
young_freq = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_freq_noERP.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_freq{i} = temp.young_freq_noERP;        % Assuming the variable is called freq_noERP
end

old_freq = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_freq_noERP.mat']);
    temp = load(file_path);             % Adjust this if the variable in the file is named differently
    old_freq{i} = temp.old_freq_noERP;        % Assuming the variable is called freq_noERP
end

% Concatenate young_freq data along the 4th dimension (rpt)
cfg = [];
young_freq = ft_appendfreq(cfg, young_freq{:}); % Concatenate all subjects
young_freq.powspctrm = permute(young_freq.powspctrm, [1, 2, 3, 4]); % Ensure subjects are in 4th dimension

% Concatenate old_freq data along the 4th dimension (rpt)
cfg = [];
old_freq = ft_appendfreq(cfg, old_freq{:}); % Concatenate all subjects
old_freq.powspctrm = permute(old_freq.powspctrm, [1, 2, 3, 4]); % Ensure subjects are in 4th dimension


%% Time-Frequency-Analysis - Power Analysis

% Step 1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'Custom_EEG1010.mat';
neighbours = ft_prepare_neighbours(cfg);

% Step 2: Configure statistics
cfg = [];
cfg.latency = [-1.5 1.5];                   % Time window for analysis
cfg.frequency = [2 100];                    % Frequency range for analysis
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
cfg.layout = 'Custom_EEG1010.mat';                                     % Specify your layout (use the same as earlier)
cfg.marker = 'on';                                          % Show electrode markers
cfg.highlight = 'on';                                       % Highlight the significant channels
cfg.highlightchannel = stat_freq.posclusterslabelmat == 1;  % Only highlight channels in the first positive cluster
cfg.colorbar = 'yes';                                       % Include colorbar for scale reference

% Plot the positive cluster = higher power in younger adults
ft_multiplotTFR(cfg, stat_freq);

% Plot MultiplotTFR for Significant Clusters
cfg = [];
cfg.layout = 'Custom_EEG1010.mat';                                                     % Use the layout specified earlier
cfg.parameter = 'stat';                                                     % Use the statistic (e.g., t-values)
cfg.colorbar = 'yes';                                                       % Include a colorbar
cfg.zlim = 'maxabs';                                                        % Adjust the color axis limits
stat_freq.mask = double(stat_freq.posclusterslabelmat == 1);                % Highlight positive cluster
cfg.maskparameter = 'mask';                                                 % Use the mask for plotting significant areas

ft_multiplotTFR(cfg, stat_freq);


% Clusterplot for Raw Power Analysis
figure
cfg=[];
cfg.layout = 'Custom_EEG1010.mat';
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
cfg.layout = 'Custom_EEG1010.mat';
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