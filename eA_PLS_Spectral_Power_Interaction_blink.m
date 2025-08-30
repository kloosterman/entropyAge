%% PLS interaction analysis for spectral power in blink data
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

%% Step 2: Load individual Spectral Power data for YA
young_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Young';

% Define abbreviations
young_ids = {'5IE98', '6AU94', '7FN98', '8RS89', '6LA93', '10PM95', '8ST94', ...
             '7SU95', '7SN94', '6FM97', '11JI91', '9JN97', '7JO97', '6KF96', ...
             '7KI87', '6MS89', '4ML96', '5SR93', '11VZ96', '10SH95'};

young_freq = cell(1, length(young_ids));
for i = 1:length(young_ids)
    file_path = fullfile(young_path, [young_ids{i}, '_freq_noERP.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    young_freq{i} = temp.young_freq_noERP;            % Assuming the variable is called mse_bin
end

% append freq data
young_freq_all = ft_appendfreq([],young_freq{:})

% Create a new copy to preserve the original
young_freq_all_demeaned = young_freq_all;

% Extract and demean powspctrm across subjects (dimension 1)
mean_pow_young = mean(young_freq_all.powspctrm, 1);  % [1 x 60 x 8]
young_freq_all_demeaned.powspctrm = young_freq_all.powspctrm - mean_pow_young;


%% Step 2: Load individual Spectral Power data for old
old_path = 'C:/Users/morit/Desktop/FoPra_Daten/Analysis_Results/Old';

% Define abbreviations
old_ids = {'8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45', '5HO51', ...
           '7MT51', '4CH63', '7CU61', '7WE49', '6RE54', '4KS63', '6JH59', ...
           '7HI40', '6CR44', '6CM48', '10SH66', '11MS53'};

old_freq = cell(1, length(old_ids));
for i = 1:length(old_ids)
    file_path = fullfile(old_path, [old_ids{i}, '_freq_noERP.mat']);
    temp = load(file_path);                 % Adjust this if the variable in the file is named differently
    old_freq{i} = temp.old_freq_noERP;            % Assuming the variable is called mse_bin
end

% append freq data
old_freq_all = ft_appendfreq([],old_freq{:})

% Create a new copy to preserve the original
old_freq_all_demeaned = old_freq_all;

% Extract and demean powspctrm across subjects (dimension 1)
mean_pow_old = mean(old_freq_all.powspctrm, 1);  % [1 x 60 x 8]
old_freq_all_demeaned.powspctrm = old_freq_all.powspctrm - mean_pow_old;


%% Step 3: Calculate group PLS analysis (Interaction)

cfg = [];
cfg.frequency = [2 100];
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
stat_freq_2group = ft_freqstatistics(cfg, young_freq_all_demeaned, old_freq_all_demeaned);


%% Step 4: Plot Interaction results

%% Step 4.1: Plot Topoplot colapsed over timescale
cfg = [];
cfg.layout = 'EEG1010';
cfg.parameter = 'stat';
cfg.interactive = 'no';
cfg.zlim = 'maxabs';
stat_freq_2group.mask = double(stat_freq_2group.stat > 3 | stat_freq_2group.stat < -3);
load colormap_jetlightgray.mat
cfg.colormap = cmap;
%cfg.maskparameter = 'mask';
cfg.colorbar = 'yes';
ft_topoplotER(cfg, stat_freq_2group);

%% Step 5.2: Plot Entropy x Timescale
%% In progress
cfg =[];
cfg.layout='EEG1010';
cfg.parameter = 'stat';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
stat_freq_2group.mask = double(stat_freq_2group.stat > 3 | stat_freq_2group.stat < -3);
load colormap_jetlightgray.mat
cfg.colormap = cmap;
cfg.maskparameter = 'mask';
ft_multiplotTFR(cfg, stat_freq_2group)

f1 = figure;
tiledlayout(2, 2); % Increase rows or columns for more subplots
cfg.clussign = 'pos';
cfg.clus2plot = 1;   cfg.integratetype = 'trapz'; % mean or trapz
stat_freq_2group.posclusterslabelmat = stat_freq_2group.mask;
ft_clusterplot3D(cfg, stat_freq_2group)

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