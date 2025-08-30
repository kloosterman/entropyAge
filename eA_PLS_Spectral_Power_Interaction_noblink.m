%% PLS interaction analysis for spectral power in noblink data
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

% append freq data
young_freq_all = ft_appendfreq([],young_freq{:})

% Create a new copy to preserve the original
young_freq_all_demeaned = young_freq_all;

% Extract and demean powspctrm across subjects (dimension 1)
mean_pow_young = mean(young_freq_all.powspctrm, 1);  % [1 x 60 x 8]
young_freq_all_demeaned.powspctrm = young_freq_all.powspctrm - mean_pow_young;


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
figure;
plot(stat_freq_2group.freq, stat_freq_2group.stat);
xlabel('Frequency (Hz)');
ylabel('Statistical Value');
title('Statistical Values over Frequency');
grid on;










%% Plot interaction results
%% Channel x Frequency

figure;
imagesc(stat_freq_2group.freq, 1:length(stat_freq_2group.label), stat_freq_2group.stat);
xlabel('Frequency (Hz)');
ylabel('Channels');
yticks(1:length(stat_freq_2group.label));
yticklabels(stat_freq_2group.label);
title('PLS Stat Values: Channels vs Frequencies');
set(gca, 'YDir', 'normal');
colormap(jet); colorbar;

%% Brainscore x Behavscore
figure;
scatter(-1*stat_freq_2group.behavscores, -1*stat_freq_2group.brainscores, 60, 'filled', ...
    'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [0 0 0]);                        % multiplied with -1 for interpretation
lsline; % regression line
xlabel('Behavioral Scores');
ylabel('Brain Scores');
r = corr(stat_freq_2group.behavscores, stat_freq_2group.brainscores, 'type', 'Spearman');
title(sprintf('Brain vs Behavior Scores (Spearman = %.2f)', r));
grid on; box on;
