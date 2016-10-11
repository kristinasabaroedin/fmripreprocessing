% This script will produce two figure
% 1) a set of figures that display the normalisation process (i.e., use to check how well normalising the EPI image to MNI space)
% 2) a motion plot that gives a visual indication of a subject's movement issues

% where spm is
% spmdir = '/media/lindenmp/WD_2TB/Dropbox/scripts/matlab/spm8/';
spmdir = '/usr/local/spm8/matlab2014a.r5236/';
% addpath(genpath(spmdir));

% where the pre-processing scripts are
% addpath(genpath('~/Dropbox/scripts/projects/OCDPG/rest_prepro/'))
addpath(genpath('/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/code/rest_prepro/'))

sublist='/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/code/sublists/trial.txt; % change to MASSIVE
fileID = fopen(sublist);
subs = textscan(fileID,'%s');
subs = subs{1};

mni_template = [spmdir,'templates/T1.nii'];

for i = 1:length(subs)
	fprintf(1,'Processing subject %s\n',subs{i})

    t1dir = ['/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/data/',subs{i},'/t1/']; 
    t1name = ['t1.nii'];  
    wm = [t1dir,'crwc2',t1name];
    csf = [t1dir,'crwc3',t1name];

    rawdir = ['/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/data/',subs{i},'/rfMRI/']; % where the unprocessed epi 4d files are
    norm_epi = ['w',epi4d];
    cd(rawdir)

    prepro_reports(mni_template,[t1dir,t1name],[t1dir,'wt1_brain.nii'],wm,csf,[rawdir,norm_epi]);

    % read in motion
    mfile = dir('rp*txt');
    mov = dlmread([rawdir,mfile(1).name]);
    mov = detrend(mov,'linear'); % detrend motion regressors

    plot_motion(mov);

    close all
end