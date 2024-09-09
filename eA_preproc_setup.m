function MEG2afc_preproc_setup()
% run from runMIBmeg_analysis

if ismac
  basepath = '/Users/kloosterman/gridmaster2012/projectdata/MEG2afc/'; %yesno or 2afc
  %     backend = 'parfor';
  backend = 'local';
  %   backend = 'qsublocal';
  compile = 'no';
else
  basepath = '/mnt/beegfs/home/kloosterman/projectdata/MEG2afc'; % on the cluster
  backend = 'slurm';
  compile = 'no';
end
stack = 1;
timreq = 240; %in minutes per run
memreq = 8000; % in MB
fun2run = @eA_preproc;

linenoise_rem = 'zapline-plus'; % bandstop zapline zapline-plus DFT
 
PREIN = fullfile(basepath, 'data');
PREOUT = fullfile(basepath, sprintf('preproc%s', linenoise_rem));

mkdir(PREOUT)
mkdir(fullfile(PREOUT, 'figures'))

overwrite = 1;

SUBJ= [1:5, 7:9, 11:21];

% study specific details
sesdirs = {'A' 'B' 'C' 'D'};
sesdir2cond = {'plac_ipsi', 'drug_ipsi', 'plac_contra', 'drug_contra'};

%make cells for each subject, to analyze in parallel
cfg = [];
cfg.PREIN = PREIN;
cfg.PREOUT = PREOUT;
cfg.linenoise_rem = linenoise_rem;
cfglist = {};

for isub = 1:length(SUBJ) 
  for ises = 1:4
    sesdir = fullfile(PREIN, ['NK' int2str(SUBJ(isub))],  sesdirs{ises}, 'meg');
    if ~exist(sesdir, 'dir'); continue; end
    cd(sesdir)
    disp(sesdir)
    meglist = dir('NK*.ds');
    for irun = 1:length(meglist)
      cfg.irun = irun;
      cfg.subjno = SUBJ(isub);
      cfg.ses = sesdirs{ises};
%       cfg.outfile = fullfile(PREOUT, sprintf('NK%d_%s_run%d.mat', SUBJ(isub), sesdir2cond{ises}, irun)); % runno appended below
      cfg.outfile = fullfile(PREOUT, sprintf('NK%d_%s_run%d_zapline.mat', SUBJ(isub), sesdir2cond{ises}, irun)); % runno appended below
      outfile_hb = fullfile(PREOUT, 'heartbeats', sprintf('NK%d_%s_run%d.mat', SUBJ(isub), sesdir2cond{ises}, irun)); % runno appended below
       
      if ~exist(outfile_hb, 'file') || overwrite
        cfglist = [cfglist cfg];
      end
    end
  end
end

% cfglist = cfglist(200)
cfglist = cfglist(randsample(length(cfglist),length(cfglist)));

fprintf('Running preproc for %d cfgs\n', length(cfglist))

if strcmp(backend, 'slurm')
  options = '-D. -c3' ; % --gres=gpu:1  
else
  options =  '-l nodes=1:ppn=3'; % torque %-q testing or gpu
end

setenv('TORQUEHOME', 'yes')
mkdir('~/qsub'); cd('~/qsub');

if strcmp(backend, 'local')
  cellfun(fun2run, cfglist)
  return
end

qsubcellfun(fun2run, cfglist, 'memreq', memreq*1e6, 'timreq', timreq*60, 'stack', stack, ...
  'StopOnError', false, 'backend', backend, 'options', options);

