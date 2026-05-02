folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

%% Turn behavior matrix into table and compute domain scores
behavfolder = '/Users/kloosterman/projectdata/EntropyAging/rsEEG-Daten/';
YA_tbl = readtable(fullfile(behavfolder, 'MA_Probanden_Neuropsychologie.xlsx'), 'Sheet', 'jüngere');
OA_tbl = readtable(fullfile(behavfolder, 'MA_Probanden_Neuropsychologie.xlsx'), 'Sheet', 'ältere');
%% Turn behavior matrix into table
load behav.mat

% todo also load excel: get age

varNames = { ...
    'd2', ...
    'VLMT_1_5', ...
    'VLMT_1', ...
    'VLMT_5', ...
    'VLMT_Interf', ...
    'VLMT_Delayed', ...
    'WMT_2', ...
    'DigitSpan', ...
    'DigitSpan_Forw', ...
    'DigitSpan_Backw', ...
    'MWT_B'};

behav_tbl = array2table(behav, 'VariableNames', varNames);

n = height(behav_tbl);

agegroup = strings(n,1);
agegroup(1:20) = "Young";
agegroup(21:end) = "Older";

behav_tbl.AgeGroup = categorical(agegroup);

%% Compute domain scores

WM_mat = [behav_tbl.DigitSpan_Forw, behav_tbl.DigitSpan_Backw];
WM_z = normalize(WM_mat);
behav_tbl.WorkingMemory = mean(WM_z, 2);

% Fluid intelligence
behav_tbl.FluidIntelligence = behav_tbl.WMT_2;

%% Z-score VLMT variables before averaging

VLMT_mat = [ ...
    behav_tbl.VLMT_1_5, ...
    behav_tbl.VLMT_1, ...
    behav_tbl.VLMT_5, ...
    behav_tbl.VLMT_Interf, ...
    behav_tbl.VLMT_Delayed];

VLMT_z = normalize(VLMT_mat);   % column-wise z-score

behav_tbl.LearningMemory = mean(VLMT_z, 2, 'omitnan');

% Attention
behav_tbl.Attention = behav_tbl.d2;

% Crystallized intelligence
behav_tbl.Crystallized = behav_tbl.MWT_B;

%% Domain table in desired order

domain_tbl = behav_tbl(:, { ...
    'AgeGroup' ...
    'WorkingMemory', ...
    'FluidIntelligence', ...
    'LearningMemory', ...
    'Attention', ...
    'Crystallized' ...
    });

%% Convert to matrix for PLS / PCA

behav_domains = table2array(domain_tbl(:,2:end));

% Reorder domains: Fluid, WM, L&M, Attention, Crystallized
ord = [2 1 3 4 5];
behav_domains = behav_domains(:,ord);

domain_labels = {'FluidIntelligence','WorkingMemory','LearningMemory','Attention','Crystallized'};

%% Optional: regimes (based on ordered domains)

% Z = normalize(behav_domains);  % z-score columns
Z = zscore(behav_domains);

% Column order now:
% 1 Fluid
% 2 WorkingMemory
% 3 LearningMemory
% 4 Attention
% 5 Crystallized

Flexible = mean(Z(:,[1 2]), 2, 'omitnan');
Balanced = mean(Z(:,[3 4]), 2, 'omitnan');
Stable   = Z(:,5);

regime_tbl = table(Flexible, Balanced, Stable);
% Regime labels
regime_labels = {'Flexible','Balanced','Stable'};

%% bar plots YA vs OA for 5 ordered z-scored domains
% Order:
% 1 Working Memory        -> Flexible
% 2 Fluid Intelligence    -> Flexible
% 3 Learning & Memory     -> Balanced
% 4 Attention             -> Balanced
% 5 Crystallized          -> Stable

% Z already exists above:
% Z = zscore(behav_domains);

% Put z-scored domains into table
domain_z_tbl = array2table(Z, 'VariableNames', ...
    {'FluidIntelligence','WorkingMemory','LearningMemory','Attention','Crystallized'});
domain_z_tbl.AgeGroup = domain_tbl.AgeGroup;

%% Remove outliers (IQR method) within each domain and group

domainNames = {'FluidIntelligence','WorkingMemory','LearningMemory','Attention','Crystallized'};
domainLabels = {'Fluid Intelligence','Working Memory','Learning & Memory','Attention','Crystallized'};

clean_tbl = domain_z_tbl;

for d = 1:length(domainNames)
    var = domainNames{d};
    
    for g = ["Young","Older"]
        idx = clean_tbl.AgeGroup == g;
        x = clean_tbl{idx, var};
        
        Q1 = prctile(x,25);
        Q3 = prctile(x,75);
        IQR = Q3 - Q1;
        
        lower = Q1 - 1.5*IQR;
        upper = Q3 + 1.5*IQR;
        
        outliers = x < lower | x > upper;
        
        % set to NaN (not delete subject globally)
        tmp = clean_tbl{idx, var};
        tmp(outliers) = NaN;
        clean_tbl{idx, var} = tmp;
    end
end

isYoung = domain_z_tbl.AgeGroup == 'Young';
isOld   = domain_z_tbl.AgeGroup == 'Older';

% Preallocate
meanY = zeros(1,5);
meanO = zeros(1,5);
semY  = zeros(1,5);
semO  = zeros(1,5);
pvals = zeros(1,5);

% Compute stats + t-tests
for i = 1:5
  % dataY = domain_z_tbl{isYoung, domainNames{i}};
  % dataO = domain_z_tbl{isOld,   domainNames{i}};
    dataY = clean_tbl{isYoung, domainNames{i}};
    dataO = clean_tbl{isOld,   domainNames{i}};

    meanY(i) = mean(dataY, 'omitnan');
    meanO(i) = mean(dataO, 'omitnan');
    
    semY(i) = std(dataY, 'omitnan') / sqrt(sum(~isnan(dataY)));
    semO(i) = std(dataO, 'omitnan') / sqrt(sum(~isnan(dataO)));
    
    [~, pvals(i)] = ttest2(dataY, dataO);
end
[p_fdr] = mafdr(pvals, 'BHFDR', true);

%% Plot
close all
% figure;  hold on
figure('Units','centimeters','Position',[5 5 8.5 7]);
 hold on
x = 1:5;
width = 0.35;

% col_Y = cmap(245,:,:);   % reddish end
% col_O = cmap(11,:,:);    % bluish end
col_Y = [1 0 0];   % reddish end
col_O =  [0 0.447 0.741];    % bluish end

% Bars
bar(x - width/2, meanY, width, ...
    'FaceColor', col_Y, 'EdgeColor', 'none', 'FaceAlpha', 0.6);

bar(x + width/2, meanO, width, ...
    'FaceColor', col_O, 'EdgeColor', 'none', 'FaceAlpha', 0.6);


jitter = 0.08; % adjust strength (0.05–0.15 works well)

for i = 1:5
    
  % Young
  x_jit_Y = (x(i)-width/2) + (rand(sum(isYoung),1)-0.5)*2*jitter;
  scatter(x_jit_Y, ...
    domain_z_tbl{isYoung, domainNames{i}}, ...
    7, col_Y, 'filled', ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.25);

  % Old
  x_jit_O = (x(i)+width/2) + (rand(sum(isOld),1)-0.5)*2*jitter;
  scatter(x_jit_O, ...
    domain_z_tbl{isOld, domainNames{i}}, ...
    7, col_O, 'filled', ...
    'MarkerFaceAlpha', 0.4, ...
    'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.25);
end

% Error bars
errorbar(x - width/2, meanY, semY, 'k', 'linestyle', 'none', 'LineWidth', 0.5);
errorbar(x + width/2, meanO, semO, 'k', 'linestyle', 'none', 'LineWidth', 0.5);

%% Significance markers
ymax_each = max([meanY + semY; meanO + semO], [], 1);
ymin_all = min([meanY - semY, meanO - semO]);
ymax_all = max([meanY + semY, meanO + semO]);
yRange = ymax_all - ymin_all;
if yRange == 0
    yRange = 1;
end
offset = 0.08 * yRange;

for i = 1:5
    % y = ymax_each(i) + offset;
    y = 2;
    
    if p_fdr(i) < 0.001
      stars = '***';
    elseif p_fdr(i) < 0.01
      stars = '**';
    elseif p_fdr(i) < 0.05
      stars = '*';
    else
      stars = 'n.s.';
    end

    plot([x(i)-width/2 x(i)+width/2], [y y], 'k', 'LineWidth', 1);
    text(x(i), y + 0.02*yRange, stars, ...
        'HorizontalAlignment', 'center', 'FontSize', 10);
end

%% Gray regime dividers
xline(2.5, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);
xline(4.5, '-', 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5);


%% Regime labels (clean version)
yl = ylim;
yr = range(yl);

regime_centers = [1.5, 3.5, 5];

% y_text = yl(1) - 0.08*yr;   % position just below x-axis
y_text = -2.6;

for i = 1:3
    text(regime_centers(i), y_text, regime_labels{i}, ...
        'HorizontalAlignment', 'center', ...
        'FontWeight', 'bold', ...
        'FontSize', 7.5);
end

% Extend lower limit to make space
ylim([yl(1) - 0.15*yr, yl(2)])

%% Formatting
xticks(x)
xticklabels(domainLabels)
xtickangle(25)
xlabel('Cognitive Domains')
ylabel('Domain score (z)')

lgd = legend({'Young','Older'}, 'Location','best');
% Control size of colored legend entries
lgd.ItemTokenSize = [10 10];   % [length height]
lgd.Position =  [  0.4766    0.8561     0.1599    0.0925];

legend boxoff

box off
set(gca, 'FontSize', 10)
ax=gca;
ax.XLim = [0.5 5.5];
ax.YLim = [-3 3];
ax.FontSize = 8;
ax.TickDir = 'out';
exportgraphics(gcf, fullfile(plotfolder,'domains_1col.pdf'), ...
    'ContentType','vector', 'BackgroundColor','white');
exportgraphics(gcf, fullfile(plotfolder,'domains_1col.png'), ...
     'BackgroundColor','white', 'Resolution', 600);


%%
% =========================
% SETTINGS
% =========================
corrtype = 'Spearman'; % 'Pearson' or 'Spearman'

% Extract only domain data (exclude AgeGroup)
X = table2array(domain_z_tbl(:,1:5));
labels = domain_z_tbl.Properties.VariableNames(1:5);

% =========================
% CORRELATION
% =========================
[R, P] = corr(X, 'Type', corrtype, 'Rows', 'pairwise');

% =========================
% PLOT
% =========================
figure;
imagesc(R)
axis square

% Color limits centered at 0
caxis([-1 1])

colormap(cmap) % or your custom map

colorbar

% Axis labels
xticks(1:5)
yticks(1:5)
xticklabels(labels)
yticklabels(labels)
xtickangle(45)

title('Correlation matrix (domains)')