folder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/stats_structs";
cd(folder)
plotfolder = "/Users/kloosterman/Library/CloudStorage/OneDrive-Personal/Documents/Entropy_in_aging/plots";

%% Turn behavior matrix into table and compute domain scores
behavfolder = '/Users/kloosterman/projectdata/EntropyAging/rsEEG-Daten/';
YA_tbl = readtable(fullfile(behavfolder, 'MA_Probanden_Neuropsychologie.xlsx'), 'Sheet', 'jüngere');
OA_tbl = readtable(fullfile(behavfolder, 'MA_Probanden_Neuropsychologie.xlsx'), 'Sheet', 'ältere');
%% Turn behavior matrix into table
load behav.mat

% todo: WM → Fluid → VLMT Learning → Attention → VLMT Interference → VLMT Delayed → Crystallized
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
domain_labels = domain_tbl.Properties.VariableNames;

%% Optional: regimes (based on ordered domains)

% Z = normalize(behav_domains);  % z-score columns
Z = zscore(behav_domains);

% Column order now:
% 1 WM
% 2 Fluid
% 3 LearningMemory
% 4 Attention
% 5 Crystallized

Flexible = mean(Z(:,[1 2]), 2, 'omitnan');
Balanced = mean(Z(:,[3 4]), 2, 'omitnan');
Stable   = Z(:,5);

regime_tbl = table(Flexible, Balanced, Stable);

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
    {'WorkingMemory','FluidIntelligence','LearningMemory','Attention','Crystallized'});
domain_z_tbl.AgeGroup = domain_tbl.AgeGroup;

%% Remove outliers (IQR method) within each domain and group

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

domainNames = {'WorkingMemory','FluidIntelligence','LearningMemory','Attention','Crystallized'};
domainLabels = {'Working Memory','Fluid Intelligence','Learning & Memory','Attention','Crystallized'};

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

% Error bars
errorbar(x - width/2, meanY, semY, 'k', 'linestyle', 'none', 'LineWidth', 1);
errorbar(x + width/2, meanO, semO, 'k', 'linestyle', 'none', 'LineWidth', 1);

% Scatter points
for i = 1:5
  scatter(repmat(x(i)-width/2, sum(isYoung), 1), ...
    domain_z_tbl{isYoung, domainNames{i}}, ...
    10, col_Y, 'filled', ...
    'MarkerFaceAlpha', 0.4, ...  'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.5);

  scatter(repmat(x(i)+width/2, sum(isOld), 1), ...
    domain_z_tbl{isOld, domainNames{i}}, ...
    10, col_O, 'filled', ...
    'MarkerFaceAlpha', 0.4, ...     'MarkerEdgeColor', 'w', ...
    'LineWidth', 0.5);
end

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
    y = ymax_each(i) + offset;
    
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

%% Regime labels
yl = ylim;
text(1.5, yl(2) + 0.06*range(yl), 'Flexible', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(3.5, yl(2) + 0.06*range(yl), 'Balanced', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
text(5.0, yl(2) + 0.06*range(yl), 'Stable', ...
    'HorizontalAlignment', 'center', 'FontWeight', 'bold');

ylim([yl(1), yl(2) + 0.12*range(yl)])

%% Formatting
xticks(x)
xticklabels(domainLabels)
xtickangle(25)
xlabel('Cognitive Domains')
ylabel('Domain score (z)')
legend({'Young','Older'}, 'Location', 'best')
legend boxoff
box off
set(gca, 'FontSize', 10)
ax=gca;
ax.XLim = [0.5 5.5];
ax.YLim = [-3 3];
ax.FontSize = 8;

exportgraphics(gcf, fullfile(plotfolder,'domains_1col.pdf'), ...
    'ContentType','vector', 'BackgroundColor','white');
exportgraphics(gcf, fullfile(plotfolder,'domains_1col.png'), ...
    'ContentType','vector', 'BackgroundColor','white');
