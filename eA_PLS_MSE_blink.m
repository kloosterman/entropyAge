%% PLS analysis MSE for both age groups separate on blink-events

% Load custom electrode + layout (created with the previous script)
load('C:/Users/morit/Desktop/FoPra_Daten/Custom_EEG1010.mat', 'elec', 'layout');

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
% Load individual Entropy data for young and old
young_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Young';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89' '4ML96', '5SR93', '11VZ96', '10SH95'};

young_mse = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_mse{i} = temp.young_mse_bin;            % Assuming the variable is called mse_bin
end

%% Step 1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'Custom_EEG1010.lay';
neighbours = ft_prepare_neighbours(cfg);

% Convert individual MSE data to freq-compatible format
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

%% Step 2: Configure statistics for PLS
cfg = [];
cfg.layout = 'Custom_EEG1010.lay'
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

%% Step 3: Compute statistics
stat_mse_young = ft_freqstatistics(cfg, young_mse_all);

% For Interpretation flip the sign for the behavscores
stat_mse_young.behavscores = -stat_mse_young.behavscores;

%% Step 4: Plot PLS results
cfg =[];
cfg.layout='Custom_EEG1010.lay';
cfg.parameter = 'stat';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
stat_mse_young.mask = double(stat_mse_young.stat > 3 | stat_mse_young.stat < -3);
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.maskparameter = 'mask';
ft_multiplotTFR(cfg, stat_mse_young)

f1 = figure;

% Upper-left subplot
nexttile;

tiledlayout(2, 2); % Increase rows or columns for more subplots
cfg.clussign = 'pos';
cfg.clus2plot = 1;   cfg.integratetype = 'trapz'; % mean or trapz
stat_mse_young.posclusterslabelmat = stat_mse_young.mask;
ft_clusterplot3D(cfg, stat_mse_young)

% Add "A" to the absolute upper-left corner of the figure
annotation('textbox', [0.05, 0.95, 0.03, 0.03], ...
    'String', 'A', 'FontSize', 18, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Bottom-left subplot
nexttile;

% Calculate Spearman correlation on original data
spearman_corr_orig_YA = corr(stat_mse_young.behavscores, stat_mse_young.brainscores, 'Type', 'Spearman');  % Spearman on the raw data

% Scatter plot with adjusted data (mean-subtracted)
s = scatter(stat_mse_young.behavscores, stat_mse_young.brainscores, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.0, 'SizeData', 40);
axis padded;
lsline;
box on;

% Adjust the x-label
xlabel('Cognitive performance', 'FontSize', 12);
xlabel_handle = xlabel('Cognitive performance', 'FontSize', 12);  % Create the x-axis label
current_position = xlabel_handle.Position;     % Get the current position of the label
xlabel_handle.Position = [current_position(1), current_position(2) - 0.1, current_position(3)]; % Move it down
ylabel('EEG Entropy at rest', 'FontSize', 12);

% Add the Spearman correlation value based on the original data in the title
title(sprintf('Spearman = %1.2f', spearman_corr_orig_YA));  % Show Spearman value from the original (untransformed) data

% Bottom-right subplot
nexttile; 

% Get the mean correlations for each cognitive measure
mean_corrs_YA = stat_mse_young.results.lvcorrs(:, 1);

% Get the upper and lower confidence intervals
ul_corrs_YA = stat_mse_young.results.boot_result.ulcorr(:, 1);
ll_corrs_YA = stat_mse_young.results.boot_result.llcorr(:, 1);

% Bar plot for correlations
b = bar(mean_corrs_YA); 
hold on;

% Set x-ticks as numbers from 1 to the number of behavioral measures
xticks(1:length(behavNames)); 
xticklabels(1:length(behavNames)); % Display numbers for behavioral measures
xlabel('Cognitive Measures','FontSize', 12);
ylabel('Correlation','FontSize', 12);
title('Brain score vs. behavior YA', 'FontSize', 12);

% Add error bars for the confidence intervals
for i = 1:length(mean_corrs_YA)
    % Calculate the height of the error bars
    err_upper_YA = ul_corrs_YA(i) - mean_corrs_YA(i);  % Distance from mean to upper bound
    err_lower_YA = mean_corrs_YA(i) - ll_corrs_YA(i);  % Distance from mean to lower bound
    
    % Use the errorbar function to plot the confidence intervals
    errorbar(i, mean_corrs_YA(i), err_lower_YA, err_upper_YA, 'k', 'LineWidth', 1.5, 'CapSize', 5);
end

set(f1, 'Position', [100, 100, 650, 500]); % 2835px wide, 2459px tall

saveas(f1, 'behavPLS_MSE_YA_Spearman', 'pdf');
saveas(f1, 'behavPLS_MSE_YA_Spearman', 'png');

%% Step 5: Calculate explained variance

eigenvalue_YA = stat_mse_young.results.s; % eigenvalues are in s
prop_variance_explained_YA = e(1)/sum(e)

%% PLS analysis MSE for OA

% Load individual Entropy data for young and old
old_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Old';

% Define abbreviations
old_ids = {'8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45', '5HO51', ...
           '7MT51', '4CH63', '7CU61', '7WE49', '6RE54', '4KS63', '6JH59', ...
           '7HI40', '6CR44', '6CM48', '10SH66', '11MS53'};

old_mse = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_mse_bin.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    old_mse{i} = temp.old_mse_bin;            % Assuming the variable is called mse_bin
end

% Convert individual MSE data to freq-compatible format
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


%% Step 2: Configure statistics for PLS
cfg = [];
cfg.layout = 'Custom_EEG1010.lay'
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 1000;                         % Number of permutation
cfg.num_boot = 1000;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 8;                            % 0 is Pearson corr, 8 is Spearman
cfg.design = behav_old;
cfg.num_cond = 1;                           % Number of conditions
cfg.num_subj_lst = size(old_mse_all.powspctrm,1); % Number of subjects per condition

%% Step 3: Compute statistics
stat_mse_old = ft_freqstatistics(cfg, old_mse_all);

% For Interpretation flip the sign for the behavscores
stat_mse_old.behavscores = -stat_mse_old.behavscores;
stat_mse_old.brainscores = -stat_mse_old.brainscores;
stat_mse_old.stat = -stat_mse_old.stat;

%% Step 4: Plot PLS results
cfg =[]; 
cfg.layout='Custom_EEG1010.lay';
cfg.parameter = 'stat';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
stat_mse_old.mask = double(stat_mse_old.stat > 3 | stat_mse_old.stat < -3);
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.maskparameter = 'mask';
ft_multiplotTFR(cfg, stat_mse_old)

f2 = figure;

% Upper-left subplot
nexttile;

tiledlayout(2,2);

% Cluster plot
cfg.clussign = 'pos';
cfg.clus2plot = 1;   
cfg.integratetype = 'trapz'; % mean or trapz
stat_mse_old.posclusterslabelmat = stat_mse_old.mask;
ft_clusterplot3D(cfg, stat_mse_old);

% Add "B" to the absolute upper-left corner of the figure
annotation('textbox', [0.05, 0.95, 0.03, 0.03], ...
    'String', 'B', 'FontSize', 14, 'FontWeight', 'bold', ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');

% Bottom-left subplot: Scatter plot with Spearman correlation and mean subtraction
nexttile;

% Calculate Spearman correlation on original data
spearman_corr_orig_OA = corr(stat_mse_old.behavscores, stat_mse_old.brainscores, 'Type', 'Spearman');  % Spearman on the raw data

% Scatter plot with adjusted data (mean-subtracted)
s = scatter(stat_mse_old.behavscores, stat_mse_old.brainscores, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0], 'LineWidth', 1.0, 'SizeData', 40);
axis padded;
lsline;
box on;

% Adjust the x-label position
xlabel('Cognitive performance');
xlabel_handle = xlabel('Cognitive performance');  % Create the x-axis label
current_position = xlabel_handle.Position;     % Get the current position of the label
xlabel_handle.Position = [current_position(1), current_position(2) - 0.1, current_position(3)]; % Move it down

ylabel('EEG Entropy at rest');
title(sprintf('Spearman = %1.2f', spearman_corr_orig_OA));  % Show Spearman value from original data

% Bottom-right subplot
nexttile; 

% Get the mean correlations for each cognitive measure
mean_corrs_OA = stat_mse_old.results.lvcorrs(:, 1);

% Get the upper and lower confidence intervals (changed because of
% changed sign)
ul_corrs_OA = stat_mse_old.results.boot_result.ulcorr(:, 1);
ll_corrs_OA = stat_mse_old.results.boot_result.llcorr(:, 1);

mean_corrs_OA = -mean_corrs_OA;
ul_corrs_OA = -ul_corrs_OA
ll_corrs_OA = -ll_corrs_OA

% Bar plot for correlations
b = bar(mean_corrs_OA); 
hold on;

% Set x-ticks as numbers from 1 to the number of behavioral measures
xticks(1:length(behavNames)); 
xticklabels(1:length(behavNames)); % Display numbers for behavioral measures
xlabel('Cognitive Measures');
ylabel('Correlation');
title('Brain score vs. behavior OA');

% Add error bars for the confidence intervals
for i = 1:length(mean_corrs_OA)
    % Calculate the height of the error bars
    err_upper_OA = ul_corrs_OA(i) - mean_corrs_OA(i);  % Distance from mean to upper bound
    err_lower_OA = mean_corrs_OA(i) - ll_corrs_OA(i);  % Distance from mean to lower bound
    
    % Use the errorbar function to plot the confidence intervals
    errorbar(i, mean_corrs_OA(i), err_lower_OA, err_upper_OA, 'k', 'LineWidth', 1.5, 'CapSize', 5);
end
set(f2, 'Position', [100, 100, 650, 500]); % 2835px wide, 2459px tall

% Save the figure in both PDF and PNG formats
saveas(f2, 'behavPLS_MSE_OA_Spearman', 'pdf');
saveas(f2, 'behavPLS_MSE_OA_Spearman', 'png');

%% Step 5: Calculate explained variance

eigenvalue_OA = stat_mse_old.results.s; % eigenvalues are in s
prop_variance_explained_OA = e(1)/sum(e)

%% t-test for brainscores

% Extract brainscores for young and old groups
young_brainscores = stat_mse_young.brainscores;
old_brainscores = stat_mse_old.brainscores;

% Perform an independent two-sample t-test
[h, p, ci, stats] = ttest2(young_brainscores, old_brainscores);

% Display results
fprintf('T-Test Results:\n');
fprintf('p-value: %.4f\n', p);
fprintf('t-statistic: %.4f\n', stats.tstat);
fprintf('Degrees of freedom: %.0f\n', stats.df);
fprintf('Confidence interval: [%.4f, %.4f]\n', ci(1), ci(2));

% Hypothesis test result
if h == 0
    disp('Fail to reject the null hypothesis: no significant difference.');
else
    disp('Reject the null hypothesis: significant difference exists.');
end