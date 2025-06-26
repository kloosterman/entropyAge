% Basic scripts to start the resting state data analysis
% This is kind of a Master-script from which to start the analysis

restoredefaultpath % reset MATLAB´s search path to the default state, this removes any previously added toolboxes to avoid conflict

addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/fieldtrip-20240916')   % set path for FieldTrip
addpath('C:/Users/morit/Desktop/Toolboxes_MATLAB/zapline-plus-main')    % set path for Zapline Toolbox
ft_defaults;                                                            % initializes Fieldtrip by setting default paths and options

% Add a path to the scripts
addpath('C:/Users/morit/Desktop/FoPra_Daten/entropyAge-main')

% Preprocessing EEG data for one subject
cfg=[];                                                                                                         % configuration structure for preprocessing script, defining variables by using cfg=[] to start definition process
cfg.dataset = 'C:/Users/morit/Desktop/FoPra_Daten/Raw_Data_Entropy_Aging/Old/9RW37/RestingState_9RW37.vhdr';  % path to the raw EEG dataset
cfg.PREOUT_FIG = 'C:/Users/morit/Desktop/FoPra_Daten';                                                          % Output path figures
cfg.PREOUT_CLEAN = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging/Young';                         % Output path clean data young
cfg.SUBJ = '9RW37';                                                                                             % Subject ID
cfglist = {};                                                                                                   % stores configuration in a cell array, which allows processing of multiple subjects if needed
cfglist{1} = cfg;                                                                                               % change the 1 to the subjects we have
cellfun(@eA_preproc_blink, cfglist)                                                                                   % applies the custom preprocessing function eA_preproc to each config in cfglist (in this case just the 7SN97)

% Perform the loop over all subjects, young and old separately

% Define the list of young subjects
%subjects = {'6AU94','7FN98','8RS89','6LA93','10PM95''8ST94','7SU95','7SN94','5PH93','6FM97', '9JN97','7JO97','7KI87','6MS89','5IE98','4ML96','5SR93','11VZ96','10SH95','10CH91'}; % Add more subjects here
% Excluded '10PM95', '11JI91', '6KF96'

% % Loop over subjects
% for i = 1:length(subjects)
%     % Define cfg for the current subject
%     cfg = [];
%     cfg.dataset = ['C:/Users/morit/Desktop/FoPra_Daten/Raw_Data_Entropy_Aging/Young/' subjects{i} '/RestingState_' subjects{i} '.vhdr'];
%     cfg.dataset_clean = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging/Young';
%     cfg.PREOUT_FIG = 'C:/Users/morit/Desktop/FoPra_Daten';
%     cfg.PREOUT_CLEAN = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging/Young';
%     cfg.SUBJ = subjects{i};
% 
%     % Call preprocessing function (e.g., eA_preproc) for each subject
%     eA_preproc_blink(cfg);
% end
% 
% % Define the list of old subjects
% subjects = {'8GH64', '8UH50', '7PE61','5MA56','5CH43','7PS49','6JU60','5BN45','5HO51','9RW37','7MT51','4CH63','7CU61','7WE49','6RE54','4KS63','6JH59','7HI40','6CR44','6CM48','10SH66','11MS53'}; % Add more subjects here
% 
% % Loop over subjects
% for i = 1:length(subjects)
%     % Define cfg for the current subject
%     cfg = [];
%     cfg.dataset = ['C:/Users/morit/Desktop/FoPra_Daten/Raw_Data_Entropy_Aging/Old/' subjects{i} '/RestingState_' subjects{i} '.vhdr'];
%     cfg.dataset_clean = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging/Old';
%     cfg.PREOUT_FIG = 'C:/Users/morit/Desktop/FoPra_Daten';
%     cfg.PREOUT_CLEAN = 'C:/Users/morit/Desktop/FoPra_Daten/Clean_Data_Entropy_Aging/Old';
%     cfg.SUBJ = subjects{i};
% 
%     % Call preprocessing function (e.g., eA_preproc) for each subject
%     eA_preproc_blink(cfg);
% end