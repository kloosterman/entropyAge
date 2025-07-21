%% PLS analysis for Entropy in noblink data
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

%% Step 2: Load individual Entropy data for YA
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

%% Step 3: PLS analysis
%% Step 3.1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'EEG1010';
neighbours = ft_prepare_neighbours(cfg);

%% Step 3.2: Configure statistics
cfg = [];
%cfg.latency = [-1.5 1.5];                   % Time window for analysis
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


%% Step 3.3: Convert individual MSE data to freq-compatible format
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

%% Step 3.4: Compute statistics
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

stat_mse_young = ft_freqstatistics(cfg, young_mse_all);

%% Step 4: Plot YA results

% General settings
cfg = [];
cfg.layout = 'EEG1010';
cfg.interactive = 'no';
cfg.zlim = 'maxabs';
load colormap_jetlightgray.mat
cfg.colormap = cmap;

%% Step 4.1 Plot Latent Variable at Timescale 100

f = find(stat_mse_young.freq == 100);
dummy = [];
dummy.label = stat_mse_young.label;
dummy.dimord = 'chan_freq_time';
dummy.freq = stat_mse_young.freq(f);
dummy.time = 0; % dummy time
dummy.powspctrm = reshape(stat_mse_young.stat(:,f), [length(stat_mse_young.label), 1, 1]);

figure;
ft_topoplotTFR(cfg, dummy);
title('YA Topomap at timescale 100');

%% Step 4.2 Plot all channels with the timescales

figure;

imagesc(stat_mse_young.freq, 1:length(stat_mse_young.label), stat_mse_young.stat);

xlabel('Timescale (Hz)');
ylabel('Channels');
yticks(1:length(stat_mse_young.label));
yticklabels(stat_mse_young.label);
set(gca, 'YDir', 'normal'); % so channel 1 is at bottom
colorbar;
title('PLS Stat Values: YA Channels vs Timescales');
colormap(jet);

%% Step 4.3: Plot behavior x brain scores
figure;
scatter(-1*stat_mse_young.behavscores, stat_mse_young.brainscores, 40, 'k', 'filled', ...
    'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.0);                                      % multiplied with -1 for the interpretation
xlabel('Behavioral Scores');
ylabel('Brain Scores');
title('Scatter plot: YA Behavior vs. Brain Scores');
grid on;
axis padded;

% Add a least-squares regression line
lsline;

% Optional: show correlation coefficient in title
r = corr(stat_mse_young.behavscores, stat_mse_young.brainscores, 'Type', 'Spearman');
title(sprintf('Behavior vs Brain Scores (Spearman = %.2f)', r));



%% Step 4.4: Plot latent variable at all timescales

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
    load colormap_jetlightgray.mat
    cfg.colormap = cmap;
    cfg.parameter = 'avg';
    cfg.comment = sprintf('Timescale %d', stat_mse_young.freq(idx));

    figure;
    ft_topoplotER(cfg, temp);
end


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

%% Step 3: PLS analysis
%% Step 3.1: Prepare neighbours
cfg = [];
cfg.method = 'triangulation';
cfg.layout = 'EEG1010';
neighbours = ft_prepare_neighbours(cfg);

%% Step 3.2: Configure statistics
cfg = [];
%cfg.latency = [-1.5 1.5];                   % Time window for analysis
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

%% Step 3.3: Convert individual MSE data to freq-compatible format
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


%% Step 3.4: Compute statistics
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
cfg.num_subj_lst = size(old_mse_all.powspctrm,1); % Number of subjects per condition

stat_mse_old = ft_freqstatistics(cfg, old_mse_all);

%% Step 4: Plot OA results

% General settings
cfg = [];
cfg.layout = 'EEG1010';
cfg.interactive = 'no';
cfg.zlim = 'maxabs';
load colormap_jetlightgray.mat
cfg.colormap = cmap;

%% Step 4.1: Plot latent variable at timescale 100
f = find(stat_mse_old.freq == 100);
dummy = [];
dummy.label = stat_mse_old.label;
dummy.dimord = 'chan_freq_time';
dummy.freq = stat_mse_old.freq(f);
dummy.time = 0; % dummy time
dummy.powspctrm = reshape(stat_mse_old.stat(:,f), [length(stat_mse_old.label), 1, 1]);

figure;
ft_topoplotTFR(cfg, dummy);
title('OA Topomap at timescale 100');

%% Step 4.2: Plot all channels with the timescale

figure;

imagesc(stat_mse_old.freq, 1:length(stat_mse_old.label), stat_mse_old.stat);

xlabel('Timescale (Hz)');
ylabel('Channels');
yticks(1:length(stat_mse_old.label));
yticklabels(stat_mse_old.label);
set(gca, 'YDir', 'normal'); % so channel 1 is at bottom
colorbar;
title('PLS Stat Values: OA Channels vs Timescales');
colormap(jet);

%% Step 4.3: Plot behavior vs. brain scores
figure;
scatter(-1*stat_mse_old.behavscores, -1*stat_mse_old.brainscores, 40, 'k', 'filled', ...
    'MarkerEdgeColor', [1 1 1], 'LineWidth', 1.0);                                      % multiplied with -1 for the interpretation
xlabel('Behavioral Scores');
ylabel('Brain Scores');
title('Scatter plot: OA Behavior vs. Brain Scores');
grid on;
axis padded;

% Add a least-squares regression line
lsline;

% Optional: show correlation coefficient in title
r = corr(stat_mse_old.behavscores, stat_mse_old.brainscores, 'Type', 'Spearman');
title(sprintf('Behavior vs Brain Scores (Spearman = %.2f)', r));

%% Step 4.4: Plot latent variable at all timescales

for idx = 1:length(stat_mse_old.freq)
    temp = [];
    temp.avg = stat_mse_old.stat(:, idx);  % data for one timescale
    temp.label = stat_mse_old.label;
    temp.dimord = 'chan_time';    % Because ft_topoplotER expects chan_time
    temp.time = 1;                % dummy time point

    cfg = [];
    cfg.layout = 'EEG1010';
    cfg.zlim = 'maxabs';
    cfg.colorbar = 'yes';
    load colormap_jetlightgray.mat
    cfg.colormap = cmap;
    cfg.parameter = 'avg';
    cfg.comment = sprintf('Timescale %d', stat_mse_old.freq(idx));

    figure;
    ft_topoplotER(cfg, temp);
end






%% Original plotting code doesn´t work
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