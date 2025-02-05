addpath('/Users/kloosterman/Documents/GitHub/entropyAge')
addpath('/Users/kloosterman/Documents/GitHub/plscmd')
addpath('/Users/kloosterman/Documents/GitHub/fieldtrip_dev')
addpath('/Users/kloosterman/Documents/GitHub/fieldtrip')
ft_defaults
load /Users/kloosterman/Documents/GitHub/plotting-tools/colormap_jetlightgray.mat
load /Users/kloosterman/Documents/GitHub/LRaudio/Acticap_64_UzL.mat; % lay comes out

datapath = '/Users/kloosterman/projectdata/EntropyAging';
cd(datapath)

behav_all = readtable('FoPra_Behavioral_Measures.xlsx');
behavNames = string(behav_all.Properties.VariableNames);
% all vars:
% behav_of_interest = ["VP"	"Sex1_m2_f"	"Age"	"Age_Group"	"Years_Education"	"HoechsterAbschlussDiplom_Master_1_Bachelor_2_Abitur_3_Real_4_Ha"	"d2KL"	"VLMTDg1_5"	"VLMTDg1"	"VLMTDg5"	"VLMTI"	"VLMTDg7"	"VLMTW"	"VLMTW_F"	"WMT_2"	"ZahlenGes"	"ZahlenVor"	"ZahlenRueck"	"MWT_B"	"MoCA"	"GDS"	"Blink_trial_all"	"Blink_trial_350"	"Blink_trial_1050"	"Blink_trial_musc_artifact"	"Blink_trial_varcut"	"Blink_Mean_Amplitude"	"Blink_Peak_Mean_Amplitude"	"Blink_Mean_Duration"];

behav_of_interest = ["Age"	"Years_Education"	"d2KL"	"VLMTDg1_5"	"VLMTDg1"	"VLMTDg5"	"VLMTI"	"VLMTDg7"	"VLMTW"	"VLMTW_F"...
                      "WMT_2"	"ZahlenGes"	"ZahlenVor"	"ZahlenRueck"	"MWT_B"];

behav_of_interest = ["d2KL"	"VLMTDg1_5"	"VLMTDg1"	"VLMTDg5"	"VLMTI"	"VLMTDg7"	"VLMTW"	"VLMTW_F"...
                      "WMT_2"	"ZahlenGes"	"ZahlenVor"	"ZahlenRueck"	"MWT_B"];

behav = behav_all(behav_all.Age_Group==1,:);
behav = behav(behav.Age_Group==1, behav_of_interest);
behav = table2array(behav);

% behav = table2array(behav(:,7:19));
% behavNames = behavNames(7:19);

% behav = behav(:,[1:6, 9:end]);
% behavNames = behavNames([1:6, 9:end]);

load('young_mse_all.mat', 'young_mse_all')

%% YA PLS analysis

cfg = [];
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 100;                         % Number of permutation
cfg.num_boot = 100;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 8;                            % 0 is Pearson corr, 8 is Spearman
cfg.num_cond = 1;                           % Number of conditions
cfg.design = behav;
cfg.num_subj_lst = size(young_mse_all.powspctrm,1); % Number of subjects per condition

% Step 3: Compute statistics
stat_mse = ft_freqstatistics(cfg, young_mse_all);

%%
cfg =[]; 
cfg.layout=lay;
cfg.parameter = 'stat';
cfg.colorbar = 'yes';
cfg.zlim = 'maxabs';
stat_mse.mask = double(stat_mse.stat > 3 | stat_mse.stat < -3);
cfg.colormap = cmap;
cfg.maskparameter = 'mask';
ft_multiplotTFR(cfg, stat_mse)

f = figure;
tiledlayout(2,2);
cfg.clussign = 'pos';
cfg.clus2plot = 1;   cfg.integratetype = 'trapz'; % mean or trapz
stat_mse.posclusterslabelmat = stat_mse.mask;
ft_clusterplot3D(cfg, stat_mse)

% plot latent behav vs latent brain?
nexttile; s = scatter(stat_mse.behavscores, stat_mse.brainscores, 'MarkerEdgeColor',[1 1 1],  'MarkerFaceColor', [0 0 0], 'LineWidth',1.0, 'SizeData', 40); axis padded; lsline; box on;
xlabel('"Cognitive rigidness"'); ylabel('EEG entropy at rest'); 
title(sprintf('r = %1.2f', corr(stat_mse.brainscores, stat_mse.behavscores)))

% bar plot corrs for each behav var
% nexttile; b=bar(corr(stat_mse.brainscores, behav)); xticklabels(behav_of_interest);
nexttile; b=bar(stat_mse.results.lvcorrs(:,1)); xticklabels(behav_of_interest);
ylabel('Correlation')
title('Brain score vs behavior')

% saveas(f, 'behavPLS_MSE', 'pdf');
% saveas(f, 'behavPLS_MSE', 'png');

%% 2 group PLS
older_mse_all = young_mse_all; % @Moritz: put the OA data here!
older_mse_all.powspctrm = older_mse_all.powspctrm(1:17,:,:,:); % pretend lower N in OA

cfg = [];
cfg.frequency = [20 100];
cfg.statistic = 'ft_statfun_pls';           % PLS statistics
cfg.num_perm = 100;                         % Number of permutation
cfg.num_boot = 100;
cfg.method = 'analytic';                    % analytic method for statistics
cfg.pls_method = 3;                         % 1 is taskPLS; 3 is behavPLS
cfg.cormode = 0;                            % 0 is Pearson corr, 8 is Spearman
cfg.num_cond = 1;                           % Number of conditions
cfg.design = [behav; behav(1:17,:)];        % append behav OA to YA
cfg.num_subj_lst = [19 17];                 % Number of subjects per condition, array!
stat_mse_2group = ft_freqstatistics(cfg, young_mse_all, older_mse_all);

% plot 2 group PLS: 1 BSR map for LV1, but two sets of bars: for YA and OA:
% Goal: see if we get same direction of bars for YA and OA
f = figure; f.Position = [680         624        1171         254];
tiledlayout(1,4);
cfg.clussign = 'pos';
cfg.clus2plot = 1;   cfg.integratetype = 'trapz'; % mean or trapz
stat_mse_2group.posclusterslabelmat = stat_mse_2group.mask;
ft_clusterplot3D(cfg, stat_mse_2group)

% plot latent behav vs latent brain
nexttile; s = scatter(stat_mse_2group.behavscores, stat_mse_2group.brainscores, 'MarkerEdgeColor',[1 1 1],  'MarkerFaceColor', [0 0 0], 'LineWidth',1.0, 'SizeData', 40); axis padded; lsline; box on;
xlabel('"Cognitive rigidness"'); ylabel('EEG entropy at rest'); 
title(sprintf('r = %1.2f', corr(stat_mse_2group.brainscores, stat_mse_2group.behavscores)))

% bar plot corrs for each behav var
nexttile; b=bar([stat_mse_2group.results.lvcorrs(1:11,1)'; stat_mse_2group.results.lvcorrs(12:end,1)']) ; xticklabels(behav_of_interest); % Assuming 11 behavioral variables!
ylabel('Correlation')
title('Brain score vs behavior YA OA')

saveas(f, 'behavPLS_MSE_2group', 'pdf');
saveas(f, 'behavPLS_MSE_2group', 'png');


%% factor analysis
behav_all = readtable('FoPra_Behavioral_Measures.xlsx');
behavNames = string(behav_all.Properties.VariableNames);

behav_of_interest = ["d2KL"	"VLMTDg1"	"VLMTDg5"	"VLMTI"	"VLMTDg7"	"VLMTW"	"VLMTW_F"...
                      "WMT_2"	"ZahlenVor"	"ZahlenRueck"	"MWT_B"];  % "VLMTDg1_5"	 "ZahlenGes"	

behav = behav_all(behav_all.Age_Group>0,:);
behav = behav(behav.Age_Group>0, behav_of_interest);
behav = varfun(@zscore, behav);

writetable(behav, 'eA_EFAdata')
behav = table2array(behav);

figure; imagesc(cov(zscore(behav))); colorbar  % = figure; imagesc(corr(behav)); colorbar

[loadings, psi, T, stats, factorScores] = factoran(behav, 2, 'scores', 'regression', 'rotate', 'promax');
biplot(loadings,'LineWidth',2,'MarkerSize',20)

YAfs = mean(factorScores(behav_all.Age_Group == 1,:))
OAfs = mean(factorScores(behav_all.Age_Group == 2,:))


