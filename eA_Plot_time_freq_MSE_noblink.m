%% Spectral Power analysis of YA and OA separately (looped) no blink

restoredefaultpath

% Add FieldTrip to MATLAB path
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916');
ft_defaults;

%% Step 1: Plotting for Older Adults
%% List of old subject IDs
old_subjects = { '8UH50', '7PE61', '5MA56', '7PS49', '6JU60', '5BN45',	'5HO51',	'7MT51',	'4CH63',	'7CU61', '7WE49',	'6RE54',	'4KS63',	'6JH59',	'7HI40',	'6CR44',	'6CM48',	'10SH66',	'11MS53'};  % Replace with actual subject IDs
data_path = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging_Controlanalysis/Old';
out_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Old';

% Load all freq_clean files
old_freq_all = {};
for i = 1:length(old_subjects)
    subj = old_subjects{i};
    file = fullfile(out_path, [subj '_data_freq_old.mat']);
    if exist(file, 'file')
        tmp = load(file);
        if isfield(tmp, 'data_freq_old') && isstruct(tmp.data_freq_old)
            old_freq_all{end+1} = tmp.data_freq_old;
        else
            warning('No valid data_freq_old in %s', file);
        end
    else
        warning('File missing: %s', file);
    end
end

%% Step 1.1: Average Spectral Power over subjects

cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';                                                  % set to 'yes' if you want to preserve subject dimension
avg_old_freq = ft_freqgrandaverage(cfg, old_freq_all{:});

% Average power over all channels (result: 1 x frequencies)
avg_old_power = squeeze(mean(avg_old_freq.powspctrm, 1, 'omitnan'));

% Get frequency values
freq_values = avg_old_freq.freq;

%% Plot 1: Average Power Spectrum 
figure;
plot(freq_values, avg_old_power, 'r-', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2)');
title('Average Power Spectrum - Older Adults');
xlim([min(freq_values), max(freq_values)]);
grid on;

%% Plot 2: Individual Averaged Spectral Power
figure;
hold on;
for i = 1:length(old_freq_all)
    subj_power = squeeze(mean(old_freq_all{i}.powspctrm, 1, 'omitnan'));  % [1 x freq]
    plot(old_freq_all{i}.freq, subj_power);
end
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2)');
title('Individual Power Spectra - Older Adults');
grid on;
legend(old_subjects, 'Location', 'northeastoutside');

% Initialize cell array
old_mse_all = cell(1, length(old_subjects));

for i = 1:length(old_subjects)
    subj = old_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % loads variable 'mse_bin'
    
    mse = tmp.mse_bin;   % shortcut
    
    % Rename 'sampen' to 'powspctrm'
    mse.powspctrm = mse.sampen;
    
    % Assign 'freq' field (FieldTrip expects this)
    mse.freq = mse.timescales;
    
    % Optionally remove 'sampen' to avoid confusion
    mse = rmfield(mse, 'sampen');
    
    % Make sure dimord is correct (should already be)
    mse.dimord = 'chan_freq_time';
    
    % Save back to the cell array
    old_mse_all{i} = mse;
end

% Initialize cell array
old_mse_all = cell(1, length(old_subjects));

for i = 1:length(old_subjects)
    subj = old_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % should contain variable 'mse_bin'
    old_mse_all{i} = tmp.mse_bin;
end

% Convert to FieldTrip-compatible freq-like structure
for i = 1:length(old_mse_all)
    old_mse_all{i}.powspctrm = old_mse_all{i}.sampen;
    old_mse_all{i}.freq = old_mse_all{i}.timescales;    % Important!
    old_mse_all{i}.dimord = 'chan_freq_time';
    old_mse_all{i} = rmfield(old_mse_all{i}, 'sampen'); % Remove original
end

%% Step 1.2: Average Entropy over subjects

cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';
avg_old_mse = ft_freqgrandaverage(cfg, old_mse_all{:});

% Result: 1 x nTimescales
avg_old_mse_value = squeeze(mean(mean(avg_old_mse.powspctrm, 1, 'omitnan'), 3, 'omitnan'));

% Extract timescales
timescales = avg_old_mse.freq;  % timescales are stored in the freq field

%% Plot 3: Average Entropy Across Timescales
figure;
plot(timescales, avg_old_mse_value, 'r-', 'LineWidth', 2);
xlabel('Timescale');
ylabel('Sample Entropy');
title('Average Entropy Across Timescales - Older Adults');
grid on;

%% Plot 4: Individual Averaged Entropy
figure;
hold on;
for i = 1:length(old_mse_all)
    subj_entropy = squeeze(mean(mean(old_mse_all{i}.powspctrm, 1, 'omitnan'), 3, 'omitnan'));  % [1 x timescales]
    plot(old_mse_all{i}.freq, subj_entropy);
end
xlabel('Timescale');
ylabel('Sample Entropy');
title('Individual Entropy Spectra - Older Adults');
grid on;
legend(old_subjects, 'Location', 'northeastoutside');


%% Step 2: Plotting for Younger Adults
%% List of young subject IDs
young_subjects = { '5IE98', '6AU94', '7FN98', '8RS89','6LA93','10PM95', '8ST94','7SU95','7SN94', '11JI91', '9JN97','7JO97', '6KF96', '7KI87' ,'6MS89', '4ML96', '5SR93','11VZ96','10SH95'};  % Replace with actual subject IDs
out_path = 'C:/Users/morit/Desktop/FoPra_Daten/Controlanalysis_Results/Young';

% Load all freq_clean files
young_freq_all = {};
for i = 1:length(young_subjects)
    subj = young_subjects{i};
    file = fullfile(out_path, [subj '_data_freq_young.mat']);
    if exist(file, 'file')
        tmp = load(file);
        if isfield(tmp, 'data_freq_young') && isstruct(tmp.data_freq_young)
            young_freq_all{end+1} = tmp.data_freq_young;
        else
            warning('No valid data_freq_young in %s', file);
        end
    else
        warning('File missing: %s', file);
    end
end

%% Step 2.1: Average Spectral Power over subjects

cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';                                                  % set to 'yes' if you want to preserve subject dimension
avg_young_freq = ft_freqgrandaverage(cfg, young_freq_all{:});

% Average power over all channels (result: 1 x frequencies)
avg_young_power = squeeze(mean(avg_young_freq.powspctrm, 1, 'omitnan'));

% Get frequency values
freq_values = avg_young_freq.freq;

%% Plot 5: Average Power Spectrum
figure;
plot(freq_values, avg_young_power, 'r-', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2)');
title('Average Power Spectrum - Younger Adults');
xlim([min(freq_values), max(freq_values)]);
grid on;

%% Plot 6: Individual Averaged Spectral Power
figure;
hold on;
for i = 1:length(young_freq_all)
    subj_power = squeeze(mean(young_freq_all{i}.powspctrm, 1, 'omitnan'));  % [1 x freq]
    plot(young_freq_all{i}.freq, subj_power);
end
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2)');
title('Individual Power Spectra - Younger Adults');
grid on;
legend(young_subjects, 'Location', 'northeastoutside');

% Initialize cell array
young_mse_all = cell(1, length(young_subjects));

for i = 1:length(young_subjects)
    subj = young_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % loads variable 'mse_bin'
    
    mse = tmp.mse_bin;   % shortcut
    
    % Rename 'sampen' to 'powspctrm'
    mse.powspctrm = mse.sampen;
    
    % Assign 'freq' field (FieldTrip expects this)
    mse.freq = mse.timescales;
    
    % Optionally remove 'sampen' to avoid confusion
    mse = rmfield(mse, 'sampen');
    
    % Make sure dimord is correct (should already be)
    mse.dimord = 'chan_freq_time';
    
    % Save back to the cell array
    young_mse_all{i} = mse;
end

% Initialize cell array
young_mse_all = cell(1, length(young_subjects));

for i = 1:length(young_subjects)
    subj = young_subjects{i};
    file = fullfile(out_path, [subj '_mse_bin.mat']);
    tmp = load(file);  % should contain variable 'mse_bin'
    young_mse_all{i} = tmp.mse_bin;
end

% Convert to FieldTrip-compatible freq-like structure
for i = 1:length(young_mse_all)
    young_mse_all{i}.powspctrm = young_mse_all{i}.sampen;
    young_mse_all{i}.freq = young_mse_all{i}.timescales;    % Important!
    young_mse_all{i}.dimord = 'chan_freq_time';
    young_mse_all{i} = rmfield(young_mse_all{i}, 'sampen'); % Remove original
end

%% Step 2.2: Average Entropy Across Timescales

cfg = [];
cfg.parameter = 'powspctrm';
cfg.keepindividual = 'no';
avg_young_mse = ft_freqgrandaverage(cfg, young_mse_all{:});

% Result: 1 x nTimescales
avg_young_mse_value = squeeze(mean(mean(avg_young_mse.powspctrm, 1, 'omitnan'), 3, 'omitnan'));

% Extract timescales
timescales = avg_young_mse.freq;  % timescales are stored in the freq field

%% Plot 7: Average Entropy Across Timescales
figure;
plot(timescales, avg_young_mse_value, 'r-', 'LineWidth', 2);
xlabel('Timescale');
ylabel('Sample Entropy');
title('Average Entropy Across Timescales - Younger Adults');
grid on;

%% Plot 8: Individual Averaged Entropy

figure;
hold on;
for i = 1:length(young_mse_all)
    subj_entropy = squeeze(mean(mean(young_mse_all{i}.powspctrm, 1, 'omitnan'), 3, 'omitnan'));  % [1 x timescales]
    plot(young_mse_all{i}.freq, subj_entropy);
end
xlabel('Timescale');
ylabel('Sample Entropy');
title('Individual Entropy Spectra - Younger Adults');
grid on;
legend(young_subjects, 'Location', 'northeastoutside');


%% Step 3: Plotting both Age Groups

%% Plot 9: Averaged Power Spectrum for Younger and Older adults
figure;
plot(freq_values, avg_old_power, 'r-', 'LineWidth', 2); hold on;
plot(freq_values, avg_young_power, 'b--', 'LineWidth', 2);
xlabel('Frequency (Hz)');
ylabel('Power (\muV^2)');
title('Average Power Spectrum: Older vs. Younger Adults');
legend({'Old', 'Young'});
grid on;

%% Plot 10: Averaged Entropy for Younger and Older Adults
figure;
plot(timescales, avg_old_mse_value, 'r-', 'LineWidth', 2); hold on;
plot(timescales, avg_young_mse_value, 'b--', 'LineWidth', 2);
xlabel('Timescale');
ylabel('Sample Entropy');
title('Entropy Across Timescales: Older vs. Younger Adults');
legend({'Old', 'Young'});
grid on;