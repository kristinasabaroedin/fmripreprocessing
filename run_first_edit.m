clear all
    %==========================================================================
    % Add paths - edit this section
    %==========================================================================

        % this section should add paths to all revevent scipts and toolboxes, which
        % include: SPM, REST, marsbar and the current set of scripts

        % where spm is
        % spmdir = '/media/lindenmp/WD_2TB/Dropbox/scripts/matlab/spm8/';
        spmdir = '/usr/local/spm8/matlab2014a.r5236/';
        % addpath(genpath(spmdir));

        % where REST is
        % addpath(genpath('~/Dropbox/scripts/matlab/REST_V1.8_130615/'))
        % directory where marsbar is. Path will be added and removed as required in
        % script below to avoid conflict with spm routines
        % marsbar_dir = '~/Dropbox/scripts/matlab/spm8/toolbox/marsbar/';

        % where the pre-processing scripts are
        % addpath(genpath('~/Dropbox/scripts/projects/OCDPG/rest_prepro/'))
        addpath(genpath('/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/code/rest_prepro/'))

        % set FSL environments 
        fsldir = '/usr/local/fsl/5.0.8/bin/'; % directory where fsl is
        setenv('FSLDIR',fsldir(1:end-4));
        setenv('FSLOUTPUTTYPE','NIFTI');

        % setup Wavelet toolbox
        % addpath('~/Dropbox/scripts/matlab/BrainWavelet')
        % setup
        
% ------------------------------------------------------------------------------
% Subject list
% ------------------------------------------------------------------------------
sublist='/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/code/sublists/trial.txt';
fileID = fopen(sublist);
subs = textscan(fileID,'%s');
subs = subs{1};

% ------------------------------------------------------------------------------
% ROIs
% ------------------------------------------------------------------------------
roidir = '/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/code/ROIspheres/';
% generate list of roi file names
roifiles = dir([roidir,'*.nii']);

% ------------------------------------------------------------------------------
% Project dir
% ------------------------------------------------------------------------------
projdir = '/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/data/';

for i = 1:length(subs)
	fprintf(1,'Processing subject %s\n',subs{i})


    % ------------------------------------------------------------------------------
    % Unzip archive
    % ------------------------------------------------------------------------------
    %subdir = [projdir,subs{i}];
    %cd(subdir)

    %jar xvf('rfMRI_Archive.zip')
    %unzip('t1_Archive.zip')

	% ------------------------------------------------------------------------------
	% Inputs
	% ------------------------------------------------------------------------------
    % Directory of processed rfMRI data
    rawdir = [projdir,subs{i},'/rfMRI/'];
    cd(rawdir)

    % Directory of T1 data
    t1dir = [projdir,subs{i},'/t1/'];
    gm = [t1dir,'crwc1t1.nii'];   

    extract_in = 'srest_prepro.nii';

    % length of time series (no. vols)
    N = 620; 
    % Repetition time of acquistion in secs
    TR = 0.754;

    %==========================================================================
    % Extract ROI time series
    %==========================================================================

    roi_ts = [];
    for i = 1:length(roifiles)
            tic;
            system([fsldir,'fslmaths ',[roidir,roifiles(i).name],' -mul ',gm,' roi_gs']);
            
            system([fsldir,'fslmeants -i ',extract_in,' -o temp.txt -m roi_gs -w']);
            
            roi_ts(:,i) = dlmread([rawdir,'temp.txt']);
            
            fprintf('Time series extracted for ROI %d of %d \n',i,length(roifiles));
            toc;
    end

    % Save out time series as .mat file
    save('DiMartino_roi_ts.mat','roi_ts')


    % ------------------------------------------------------------------------------
    % Save relevant variables
    % ------------------------------------------------------------------------------

    % Left hemi
    R = roi_ts(:,1:6);
    save spm_regs_L R

    % Right hemi
    R = roi_ts(:,7:12);
    save spm_regs_R R

	% ------------------------------------------------------------------------------
	% Run first level
	% ------------------------------------------------------------------------------
    

    cd(rawdir);

    % Split EPI into 3D files
    spm_file_split(extract_in)

    % estimate first level for left hemisphere
    mkdir('FirstLevel_L_DiMartino');
    movefile('spm_regs_L.mat','FirstLevel_L_DiMartino')
    
    first_level_3D([rawdir,'FirstLevel_L_DiMartino/'],'scans',TR,[rawdir,extract_in],N,[rawdir,'FirstLevel_L_DiMartino/spm_regs_L.mat'],[t1dir,'wt1_brain.nii']);
    % run contrasts
    first_level_contrasts(0,[rawdir,'FirstLevel_L_DiMartino/SPM.mat'])

   
    % estimate first level for right hemisphere
    mkdir('FirstLevel_R_DiMartino');
    movefile('spm_regs_R.mat','FirstLevel_R_DiMartino')
    
    first_level_3D([rawdir,'FirstLevel_R_DiMartino/'],'scans',TR,[rawdir,extract_in],N,[rawdir,'FirstLevel_R_DiMartino/spm_regs_R.mat'],[t1dir,'wt1_brain.nii']);
    % run contrasts
    first_level_contrasts(0,[rawdir,'FirstLevel_R_DiMartino/SPM.mat'])

    system('rm -f srest_prepro_*.nii')

    fprintf('First level analysis done \n');

    

    % ------------------------------------------------------------------------------
    % rezip
    % ------------------------------------------------------------------------------
    % % Zip up rfMRI outputs

    % % Change in subject's directory (i.e., one above rfMRI directory)
    % cd(rawdir)
    % cd ../

    % % Make output directory for first level (we want to retain an unzipped copy of this)
    % mkdir('rfMRI_FirstLevel_DiMartino')
    % % Copy First level subdirs
    % copyfile([rawdir,'FirstLevel_L_DiMartino'],'rfMRI_FirstLevel_DiMartino/FirstLevel_L')
    % copyfile([rawdir,'FirstLevel_R_DiMartino'],'rfMRI_FirstLevel_DiMartino/FirstLevel_R')

    % Make output directory for preprocessed nii file (we want to retain an unzipped copy of this)
    % mkdir('rfMRI_preprocessed')
    % % Copy preprocessed nifti file
    % copyfile(rawdir,'srest_prepro.nii')
    

    % % Zip rfMRI and t1 dir
    % zip('rfMRI_Archive.zip','rfMRI')
    % zip('t1_Archive.zip','t1')

    % % Delete rfMRI and t1 dir
    % rmdir('rfMRI','s')
    % rmdir('t1','s')    

end


