% PLS
restoredefaultpath;

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

% Add path to the behavioral measures file
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

% Load individual Entropy data for YA
young_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Young';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89', '4ML96', '5SR93', '11VZ96', '10SH95'};

young_mse = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_mse{i} = temp.mse_bin;            % Assuming the variable is called mse_bin
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
cfg.num_subj_lst = size(young_mse_all.powspctrm,1); % Number of subjects per condition

%% Step 5: Compute statistics
stat_mse_young = ft_freqstatistics(cfg, young_mse_all);

%% Plot
cfg =[]; 
cfg.layout='EEG1010';
cfg.parameter = 'stat';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
stat_mse_young.mask = double(stat_mse_young.stat > 3 | stat_mse_young.stat < -3);
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.maskparameter = 'mask';
ft_multiplotTFR(cfg, stat_mse_young);

f = figure;
tiledlayout(2,2);
cfg.clussign = 'pos';
cfg.clus2plot = 1;   cfg.integratetype = 'trapz'; % mean or trapz
stat_mse_young.posclusterslabelmat = stat_mse_young.mask;
ft_clusterplot3D(cfg, stat_mse_young)

% plot latent behav vs latent brain?
nexttile; s = scatter(stat_mse_young.behavscores, stat_mse_young.brainscores, 'MarkerEdgeColor',[1 1 1],  'MarkerFaceColor', [0 0 0], 'LineWidth',1.0, 'SizeData', 40); axis padded; lsline; box on;
xlabel('"Cognitive rigidness"'); ylabel('EEG entropy at rest'); 
title(sprintf('r = %1.2f', corr(stat_mse_young.brainscores, stat_mse_young.behavscores)))

% bar plot corrs for each behav var
% nexttile; b=bar(corr(stat_mse.brainscores, behav)); xticklabels(behav_of_interest);
nexttile; b=bar(stat_mse_young.results.lvcorrs(:,1)); xticklabels(behavNames);
ylabel('Correlation')
title('Brain score vs behavior')

% saveas(f, 'behavPLS_MSE', 'pdf');
% saveas(f, 'behavPLS_MSE', 'png');



%% All timescale in topo separately

for idx = 1:length(stat_mse_young.freq)
    temp = [];
    temp.avg = stat_mse_young.stat(:, idx);  % data for one timescale
    temp.label = stat_mse_young.label;
    temp.dimord = 'chan_time';    % Because ft_topoplotER expects chan_time
    temp.time = 1;                % dummy time point

    cfg = [];
    cfg.layout = 'EEG1010';
    cfg.zlim = 'maxabs';
    cfg.colorbar = 'yes';
    cfg.parameter = 'avg';
    cfg.comment = sprintf('Timescale %d', stat_mse_young.freq(idx));

    figure;
    ft_topoplotER(cfg, temp);
end








%%
cfg =[];
cfg.layout='EEG1010';
cfg.parameter = 'stat';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
stat_mse_young.mask = double(stat_mse_young.stat > 3 | stat_mse_young.stat < -3);
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.maskparameter = 'mask';
ft_multiplotTFR(cfg, stat_mse_young)

f1 = figure;
tiledlayout(2, 2); % Increase rows or columns for more subplots
cfg.clussign = 'pos';
cfg.clus2plot = 1;   cfg.integratetype = 'trapz'; % mean or trapz
stat_mse_young.posclusterslabelmat = stat_mse_young.mask;
ft_clusterplot3D(cfg, stat_mse_young)

% Bottom-left subplot
nexttile;

% Original data (without mean subtraction)
x_data_orig = stat_mse_young.behavscores;
y_data_orig = stat_mse_young.brainscores;

% Calculate Spearman correlation on original data
spearman_corr_orig = corr(x_data_orig, y_data_orig, 'Type', 'Spearman');  % Spearman on the raw data

% Subtract the mean from both x and y for visualization purposes
x_data = x_data_orig - mean(x_data_orig);
y_data = y_data_orig - mean(y_data_orig);

% Scatter plot with adjusted data (mean-subtracted)
s = scatter(x_data, y_data, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.0, 'SizeData', 40);
axis padded;
lsline;
box on;

% Adjust the x-label
xlabel('Cognitive rigidness');
xlabel_handle = xlabel('Cognitive rigidness');  % Create the x-axis label
current_position = xlabel_handle.Position;     % Get the current position of the label
xlabel_handle.Position = [current_position(1), current_position(2) - 0.1, current_position(3)]; % Move it down

ylabel('EEG Entropy at rest');

% Add the Spearman correlation value based on the original data in the title
title(sprintf('Spearman = %1.2f', spearman_corr_orig));  % Show Spearman value from the original (untransformed) data

% Get the position of the bottom-left subplot (in normalized figure units)
ax = gca;  % Get current axis
subplot_position = ax.Position;  % Position is in [left, bottom, width, height] normalized to [0, 1]

% Right arrow pointing to the right
annotation('textarrow', ...
    [subplot_position(1) + 0.15, subplot_position(1) + 0.225], ... % x position (left to right)
    [subplot_position(2) - 0.035, subplot_position(2) - 0.035], ... % y position (slightly below the plot)
    'Color', 'black', 'LineWidth', 1.5);

% Left arrow pointing to the left
annotation('textarrow', ...
    [subplot_position(1) + 0.15, subplot_position(1) + 0.075], ... % x position (right to left)
    [subplot_position(2) - 0.035, subplot_position(2) - 0.035], ... % y position (slightly below the plot)
    'Color', 'black', 'LineWidth', 1.5);

%Bottom-left subplot
% nexttile; s = scatter(stat_mse_young.behavscores, stat_mse_young.brainscores, 'MarkerEdgeColor',[1 1 1],  'MarkerFaceColor', [0 0 0], 'LineWidth',1.0, 'SizeData', 40); axis padded; lsline; box on;
% xlabel('Cognitive rigidness'); ylabel('EEG Entropy at rest')
% title(sprintf('Spearman = %1.2f', corr(stat_mse_young.brainscores, stat_mse_young.behavscores)))

% Bottom-right subplot
nexttile; 

% Get the mean correlations for each cognitive measure
mean_corrs = stat_mse_young.results.lvcorrs(:, 1);

% Get the upper and lower confidence intervals
ul_corrs = stat_mse_young.results.boot_result.ulcorr(:, 1);
ll_corrs = stat_mse_young.results.boot_result.llcorr(:, 1);

% Bar plot for correlations
b = bar(mean_corrs); 
hold on;

% Set x-ticks as numbers from 1 to the number of behavioral measures
xticks(1:length(behavNames)); 
xticklabels(1:length(behavNames)); % Display numbers for behavioral measures
xlabel('Cognitive Measures');
ylabel('Correlation');
title('Brain score vs. behavior YA');

% Add error bars for the confidence intervals
for i = 1:length(mean_corrs)
    % Calculate the height of the error bars
    err_upper = ul_corrs(i) - mean_corrs(i);  % Distance from mean to upper bound
    err_lower = mean_corrs(i) - ll_corrs(i);  % Distance from mean to lower bound
    
    % Use the errorbar function to plot the confidence intervals
    errorbar(i, mean_corrs(i), err_lower, err_upper, 'k', 'LineWidth', 1.5, 'CapSize', 5);
end

saveas(f1, 'behavPLS_MSE_YA_Spearman', 'pdf');
saveas(f1, 'behavPLS_MSE_YA_Spearman', 'png');