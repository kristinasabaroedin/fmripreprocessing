function [] = run_all_first(subject)

    % This script will perform integrated pre-processing and ROI time series
    % extraction for a single subject. It can easily be converted into a
    % function to loop over subjects. It can be divided into two
    % streams--one for T1 and the other for EPI processing.
    %
    % Only some basic inputs need to be defined at the beginning of the script.
    %
    % The nuisance regression is largely based on a variant of the CompCor
    % procedure first described by Behzadi et al. (2007) NeuroImage, which was
    % shown to provide noise correction to fMRI data comparable to RETROICOR. A
    % subsequent variant, described by Muschelli et al. (2014) NeuroImage, was
    % suggested to also offer good correction for head motion, with no
    % additional gains provided by scrubbing.
    %
    % NB: The latest version of FSL, which allows the '-w' option in fslmeants, 
    % should be installed
    %
    % 1 - Segment T1 (output native and normalized tissue masks)
    % 2 - Slice-timing correction (optional)
    % 3 - Spatially normalize T1 to MNI template
    % 4 - EPI Realignment
    % 5 - Co-registration of realigned images to T1 (without resampling) 
    % 6 - Application of T1-spatial normalization parameters to coregistered
    %     epi (Note: bounding box has been changed to output data to 91x109x91)
    % 7 - Generate WM and CSF masks and optional mask erosion for CompCor
    % 8 - Linear detrending of realigned EPI time series (uses REST)
    % 9 - Voxel-wise time series extraction (uses fslmeants) from nuisance
    %     masks
    % 10 - CompCor analysis of wm/csf voxel time course to generate noise
    %      regressors
    % 11 - Nuisance regression of realigned detrended EPI time courses against CompCor 
    %      regressors generated in step 10 + 6 head motion parameters and their
    %      derivatives (uses fsl_regfilt)
    % 12 - Despike with Patel's Brain WaveLet toolbox (optional - not working yet...)
    % 13 - Bandpass filtering of cleaned data (uses REST)
    % 14 - Spatially smooth filtered EPIs (optional)
    % 15 - Extract ROI time series (uses fslmeants). Time series will be
    %      weighted by GM probability
    % 16 - Print diagnostic reports
    %
    % *************************************************************************
    %       NB: You must first change line 73 in spm_defaults.m to:
    %                   defaults.mask.thresh = -Inf
    %                 for first level analysis to run
    % *************************************************************************                   
    % 
    %
    % Linden Parkes, June 2015.
    %
    %==========================================================================

    % clear all;
    % subject = '1008.2.48.9';

    cnt = 1;

    %==========================================================================
    % Add paths - edit this section
    %==========================================================================

        % this section should add paths to all revevent scipts and toolboxes, which
        % include: SPM, REST, marsbar and the current set of scripts

        % where spm is
        spmdir = '/media/lindenmp/WD_2TB/Dropbox/scripts/matlab/spm8/';
        addpath(genpath(spmdir));

        % where REST is
        addpath(genpath('~/Dropbox/scripts/matlab/REST_V1.8_130615/'))
        % directory where marsbar is. Path will be added and removed as required in
        % script below to avoid conflict with spm routines
        % marsbar_dir = '~/Dropbox/scripts/matlab/spm8/toolbox/marsbar/';

        % where the pre-processing scripts are
        addpath(genpath('~/Dropbox/scripts/projects/OCDPG/rest_prepro/'))

        % set FSL environments 
        fsldir = '/usr/share/fsl/5.0/bin/'; % directory where fsl is
        setenv('FSLDIR',fsldir(1:end-4));
        setenv('FSLOUTPUTTYPE','NIFTI');

        % setup Wavelet toolbox
        addpath('~/Dropbox/scripts/matlab/BrainWavelet')
        setup

    %==========================================================================
    % Basic inputs - edit this section 
    %==========================================================================

        % directory where raw .nii files are. 
        rawdir = ['/media/lindenmp/SG8_4TB/Research_Projects/OCDPG/data/',subject,'/rfMRI/']; % where the unprocessed epi 4d files are
        % file name of EPI 4d file
        epifile = dir([rawdir,'epi.nii']);
        epi4d = epifile(1).name;

        % length of time series (no. vols)
        N = 189; 
        % Repetition time of acquistion in secs
        TR = 2.5;
        % Desired voxel dimension (in mm) of analysis after spatial normalization
        voxdim = 2;

        % Do you want to run slice-timing correction? 1 = yes; 0 = no.
        slicetime = 1;
        % Number of slices in epi volumes. 
        % Set to [] if slicetime = 0.
        nslices = 44;
        % Vector defining acquisition order of EPI slices (necessary for slice-timing correction.
        % See help of slicetime_epis.m for guidance on how to define)
        % Set to [] if slicetime = 0.
        % order = [1:1:nslices]; % ascending
        order = [1:2:nslices-1,2:2:nslices]; % interleaved
        % Reference slice for slice timing acquisition. See help of slicetime_epis.m 
        % for guidance on how to define. 
        % Set to [] if slicetime = 0.
        refslice = nslices-1;

        % Low-pass cut-off for bandpass filter in Hz (e.g., .08) 
        LoCut = 0.08;
        % Hi-pass cut-off for bandpass filter in Hz (e.g., .008)
        HiCut = 0.008;



        % Directory where the t1 is
        t1dir = ['/media/lindenmp/SG8_4TB/Research_Projects/OCDPG/data/',subject,'/t1/']; 
        % name of t1 file.
        t1name = ['t1.nii'];  

        % the path and filename of the template in MNI space to which everything
        % will be normalized
        mni_template = [spmdir,'templates/T1.nii'];
        %mni_template = [fsldir(1:end-4),'/data/standard/MNI152_T1_2mm.nii'];

        % directory where roi .mat files generated by marsbar are. Note - there
        % should be not other .mat files in this directory. Only the roi .mat files
        roidir = '/media/lindenmp/SG8_4TB/Research_Projects/OCDPG/TriStri/';

        % generate list of roi file names
        roifiles = dir([roidir,'*.nii']);
        % Remove the all-in-one ROI files, we only want the separate binary masks.
        roifiles(strcmp({roifiles.name},'Lk3.nii')) = [];
        roifiles(strcmp({roifiles.name},'Rk3.nii')) = [];

        % Compute 1st-level stats? 1 = yes; 0 = no.
        run_stats = 0;

        % OPTIONAL STEPS

        % Option to erode wm and csf noise masks by 1 voxel. 0 = no; 1 = yes.
        erosion = 0;

        % Wavelet despiking (Patel 2014 Neuroimage)
        waveDespike = 0;

        % Scalar value indicating spatial smoothing kernal size in mm. 
        % To skip smoothing, set to kernel = [];
        kernel = 8;

    %==========================================================================
    % Preprocess T1
    % In this step, we take the subject's T1 image and segment it into the three
    % basic tissue types using SPM8
    % The tissue types are:
    % 1) grey matter
    % 2) white matter
    % 3) cerebrospinal fluid
    % The results are probabilistic brain images where higher values in voxels
    % represent a greater chance that a given voxel is a given tissue type
    %==========================================================================

        tic;

        cd(t1dir);

        % Tissue segment T1 with SPM
        % segment_t1([t1dir,t1name],spmdir);
         
        % Define outputs
        gm = [t1dir,'crwc1',t1name];
        wm = [t1dir,'crwc2',t1name];
        csf = [t1dir,'crwc3',t1name];
         
        fprintf('T1 segmented \n');

        times(cnt) = toc;
        cnt = cnt+1;

    %==========================================================================
    % Slice timing correction, realignment, and normalization
    % In this step, we perform some standard steps for pre-processing functional
    % data using SPM8
    %==========================================================================

        cd(rawdir);

        tic;
        % Run optional slice-timing correction and set inputs and outputs of
        % realignment accordingly
        if slicetime == 1 
            
            % slice-timing correction
            slicetime_epis([rawdir,epi4d], nslices, TR, order, refslice, N); 

            % I/O names
            realign_in = ['a',epi4d];
            realign_out = ['ra',epi4d];
            norm_epi = ['wa',epi4d];
            
        else
            
            % I/O names if slice-time correction not performed
            realign_in = epi4d ;
            realign_out = ['r',epi4d];
            norm_epi = ['w',epi4d];
            
        end

        % realignment
        spatial_prepro_4d([rawdir,realign_in],[t1dir,t1name],mni_template,voxdim,N);  


        fprintf('Slice-timing correction, realignment and EPI to T1 co-registration done \n');
        times(cnt) = toc;
        cnt = cnt+1;

    %==========================================================================
    % Generate tissue masks
    % In this step we generate two mask that we use later on as nuisance regressors.
    % These are masks of the white matter and csf
    % They are then used with CompCor
    %==========================================================================

        tic;

        cd(t1dir);

        % Threshold GM and WM tissue maps at .10 and csf at .50
        system([fsldir,'fslmaths ',gm,' -thr .10 gm_thresh']);
        system([fsldir,'fslmaths ',wm,' -thr .10 wm_thresh']);
        system([fsldir,'fslmaths ',csf,' -thr .50 csf_thresh']);
         
        % Add GM and WM to create a brain mask
        system([fsldir,'fslmaths gm_thresh -add wm_thresh -add csf_thresh -bin ',t1dir,'t1_segmask']);
        system([fsldir,'fslmaths ',t1dir,'w',t1name,' -mul t1_segmask wt1_brain']);   % multiply by t1 to have non-binary version

        % threshold images
        system([fsldir,'fslmaths ',gm,' -thr .01 gm01']);  % retain anything with >=1% probability of being gm
        system([fsldir,'fslmaths ',wm,' -thr .99 wm99']); % retain only voxels with >=99% probability of wm
        system([fsldir,'fslmaths ',csf,' -thr .99 csf99']); % retain only voxels with >=99% probability of csf

        % binarize thresholded images
        system([fsldir,'fslmaths gm01 -bin gm01_bin'])
        system([fsldir,'fslmaths wm99 -bin wm99_bin'])
        system([fsldir,'fslmaths csf99 -bin csf99_bin'])

        % remove overlap between gm and white and csf masks
        system([fsldir,'fslmaths gm01_bin -mul -1 -add 1 gm01_inv']) ;
        system([fsldir,'fslmaths wm99_bin -mul gm01_inv wm_final']) ;
        system([fsldir,'fslmaths csf99_bin -mul gm01_inv csf_final']) ;

        % multiply to get back probability values
        system([fsldir,'fslmaths wm_final -mul ',wm,' wm_final']) ;
        system([fsldir,'fslmaths csf_final -mul ',csf,' csf_final']) ;

        % erode wm and csf masks if chosen
        if erosion == 1
            
            system([fsldir,'fslmaths wm_final.nii -ero ',t1dir,'wm_final_erode']);
            system([fsldir,'fslmaths csf_final.nii -ero ',t1dir,'csf_final_erode']);
            system('gunzip -f *erode*');
            
            wmfinal = 'wm_final_erode.nii';
            csffinal = 'csf_final_erode.nii';
            
        else
            
            wmfinal = 'wm_final.nii';
            csffinal = 'csf_final.nii';
            
        end

        fprintf('Nuisance masks generated \n');
        times(cnt) = toc;
        cnt = cnt+1;

    % ------------------------------------------------------------------------------
    % Mask out non-brain tissue
    % ------------------------------------------------------------------------------

        fprintf(1, 'Extracting brain..\n');
        % mask filtered volume to retain only brain
        system([fsldir,'fslmaths ',rawdir,norm_epi,' -mas ',t1dir,...
            'wt1_brain ',rawdir,'brain_',norm_epi]);

    %==========================================================================
    % Detrend epis
    %==========================================================================

        DetrendIn = ['brain_',norm_epi];

        tic;

        cd(rawdir)

        system(['mkdir temp']); % create separate directory for REST detrend function
        dtdir=[rawdir,'temp/']; 
        system(['mv ',rawdir,DetrendIn,' ',dtdir]); % move 4d file to directory

        rest_detrend(dtdir, '_detrend') % detrend epis using rest
        % This write out a file called 'detrend_4DVolume'

        system('mv temp/* .'); % move files back to rawdir
        system('mv temp_detrend/* .'); % move files back to rawdir

        system('rmdir temp temp_detrend'); % delete directories

        fprintf('Detrending complete \n');
        times(cnt) = toc;
        cnt = cnt+1;


    % ------------------------------------------------------------------------------
    % Patel's Wavelet despiking
    % ------------------------------------------------------------------------------

        if waveDespike == 1
            cd(rawdir)
            
            WaveletDespike('detrend_4DVolume.nii','waveDespike','LimitRAM',10);

            system('gunzip -f waveDespike*.gz');
        end


    %==========================================================================
    % Extract nuisance time courses and run CompCor
    %==========================================================================

        if waveDespike
            CompCorIn = 'waveDespike_wds.nii';
        else
            CompCorIn = 'detrend_4DVolume.nii';
        end
            

        tic;

        cd(rawdir);
         
        system([fsldir,'fslmeants -i ',CompCorIn,' -o wm_ts.txt -m ',t1dir,wmfinal,' --showall']); 
        system([fsldir,'fslmeants -i ',CompCorIn,' -o csf_ts.txt -m ',t1dir,csffinal,' --showall']); 

        % Read in wm time courses
        TCwm = dlmread('wm_ts.txt');
        TCwm(1:3,:) = []; % delete first 3 rows, which list voxel coords

        % run wm compcor
        vecs_wm = compcor(TCwm,'retain',[],5);
        clear TCwm % clear to save memory

        % Read in csf time courses
        TCcsf = dlmread('csf_ts.txt');
        TCcsf(1:3,:) = []; % delete first 3 rows, which list voxel coords

        % run csf compcor
        vecs_csf = compcor(TCcsf,'retain',[],5);
        clear TCcsf % clear to save memory

        fprintf('CompCor complete \n'); 
        times(cnt) = toc;
        cnt = cnt+1;

    %==========================================================================
    % Clean data: fsl_regfilt
    %==========================================================================

        tic;

        cd(rawdir)

        % read in motion
        mfile = dir('rp*txt');
        mov = dlmread([rawdir,mfile(1).name]);
        mov = detrend(mov,'linear'); % detrend motion regressors
        mov_diff = [zeros(1,size(mov,2)); diff(mov)]; % compute differentials

        % Combine all motion regressors and detrend
        all_mov = [mov mov_diff];

        % combine all noise
        all = [all_mov vecs_wm vecs_csf];

        % write out noise signals as text file
        dlmwrite('noise_signals.txt',all,'delimiter','\t','precision','%.6f');

        % clean data with fsl_regfilt
        % Linden: x is a variable that stores the -f flag input for fsl_regfilt.
        x = regexprep(num2str(1:size(all,2)),' ',',');
        x = regexprep(x,',,,',',');
        x = regexprep(x,',,',',');
        x = ['"',x,'"'];

        system([fsldir,'fsl_regfilt -i ',CompCorIn,' -o detrend_clean -d noise_signals.txt -f ',x]);
        % This write out a file called 'detrend_clean'

        clear x

        system('gunzip -rf *detrend_clean*');

        fprintf('Nuisance regression complete \n'); 
        times(cnt) = toc;
        cnt = cnt+1;

    %==========================================================================
    % Bandpass filter with REST
    %==========================================================================

        tic;

        % bandpass filter with rest toolbox
        cd(rawdir)

        system('mkdir temp'); % creater input directory for REST function
        cleandir=[rawdir,'temp/'];
        system(['mv ',rawdir,'detrend_clean.nii ',cleandir]);

        rest_bandpass(cleandir,TR,LoCut,HiCut,'No',0,1) % bandpass epis using rest
        % This write out a file called 'Filtered_4DVolume'

        system('mv temp/* .');         % move files back into rawdir
        system('mv temp_filtered/* .');
        system('rmdir temp temp_filtered'); % delete directories

        % Rename intuitively!!
        movefile('Filtered_4DVolume.nii','rest_prepro.nii')

        % % mask filtered volume to retain only brain
        % system([fsldir,'fslmaths ',rawdir,'Filtered_4DVolume -mas ',t1dir,...
        %     'wt1_brain ',rawdir,'Filtered_4DVolume_brain']);

        fprintf('Bandpass filtering complete \n'); 
        times(cnt) = toc;
        cnt = cnt+1;

    %==========================================================================
    % Spatially smooth the data (optional)
    %==========================================================================

        if ~isempty(kernel)
            tic;
            
            smooth_epi([rawdir,'rest_prepro.nii'],kernel,N);

            fprintf('Spatial smoothing complete \n');
            times(cnt) = toc;
            cnt = cnt+1;
        end

    %==========================================================================
    % Extract ROI time series
    %==========================================================================

        tic;

        cd(rawdir)

        % defin input name depending on whether smoothing was performed
        if ~isempty(kernel)
            extract_in = 'srest_prepro.nii';
        else
            extract_in = 'rest_prepro.nii';
        end

        roi_ts = [];

        %generate list of rois
        for i = 1:length(roifiles)
            
            tic;
            system([fsldir,'fslmaths ',[roidir,roifiles(i).name],' -mul ',gm,' roi_gs']);
            
            system([fsldir,'fslmeants -i ',extract_in,' -o temp.txt -m roi_gs -w']);
            
            roi_ts(:,i) = dlmread([rawdir,'temp.txt']);
            
            fprintf('Time series extracted for ROI %d of %d \n',i,length(roifiles)); toc;
            
        end

        system(['rm ',rawdir,'temp.txt']);
        system(['rm ',rawdir,'roi_gs*']);

        fprintf('ROI time series extracted \n'); 
        times(cnt) = toc;
        cnt = cnt+1;

    %==========================================================================
    % Print reports
    %==========================================================================

        tic;
        
        plot_motion(mov);

        prepro_reports(mni_template,[t1dir,t1name],[t1dir,'wt1_brain.nii'],wm,csf,[rawdir,norm_epi]);

        times(cnt) = toc;
        cnt = cnt+1;

        close all

    %==========================================================================
    % Save relevant variables
    %==========================================================================

        tic;
         
        cd(rawdir);

        save roi_ts roi_ts
        save noise_ts vecs_wm vecs_csf *mov* all

        % this step has been modified to assume that you have an even number of ROIs and then they are split into two halves.
        % the first half are all from the left hemisphere and the second half are all from the right hemisphere.
        % in my case, this translates to 3 Left hemi ROIs and 3 right hemi ROIs - Linden 10/6/15
        R = roi_ts(:,1:length(roifiles)/2);
        save spm_regs_L R
        clear R

        R = roi_ts(:,length(roifiles)/2+1:end);
        save spm_regs_R R
        clear R
        % Original code - i.e., take all the ROIs
        % R = roi_ts;
        % save spm_regs R

        fprintf(' Processing complete! \n Total time taken: %d mins \n Number of noise regressors: %d \n',...
            sum(times)/60, size(all,2));

        times(cnt) = toc;
        cnt = cnt+1;

    %==========================================================================
    % Run first level specification and estimation
    %==========================================================================

        if run_stats
            
            tic;
            cd(rawdir);

            % estimate first level for left hemisphere
            mkdir('FirstLevel_L');
            movefile('spm_regs_L.mat','FirstLevel_L')
            first_level([rawdir,'FirstLevel_L/'],'scans',TR,[rawdir,extract_in],N,[rawdir,'FirstLevel_L/spm_regs_L.mat'],[t1dir,'wt1_brain.nii']);
            % run contrasts
            first_level_contrasts(1,[rawdir,'FirstLevel_L/SPM.mat'])

            % estimate first level for right hemisphere
            mkdir('FirstLevel_R');
            movefile('spm_regs_R.mat','FirstLevel_R')
            first_level([rawdir,'FirstLevel_R/'],'scans',TR,[rawdir,extract_in],N,[rawdir,'FirstLevel_R/spm_regs_R.mat'],[t1dir,'wt1_brain.nii']);
            % run contrasts
            first_level_contrasts(1,[rawdir,'FirstLevel_R/SPM.mat'])


            fprintf('First level analysis done \n'); toc;
            times(cnt) = toc;
            cnt = cnt+1;
            
        end

    save times times
end