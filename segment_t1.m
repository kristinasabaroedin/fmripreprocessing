function [] = segment_t1(t1file,spmdir)

% segment_t1(t1file,spmdir)
%
% This function will perform tissue segmentation of a t1 with SPM using 
% the SPM8 New Segment routine. It will also resample the segmented tissue
% masks to MNI dimensions.
%
% -------
% INPUTS:
% -------
%
% t1file   - a string containing the full path and filename of the T1. 
%            e.g., '/path/to/t1/t1.nii'
%
% spmdir    - string contating the path to where spm is installed; e.g.,
%           '/usr/local/spm8'.
%
% -------
% OUTPUTS:
% -------
%
% segmented images will be written with the prefix ‘c’. e.g., ‘c1’ refers to gm,
% ‘c2’ to wm, ‘c3’ to csf.
%
% normalised, unmodulated images will be written with the prefix ‘wc’; modulated 
% images with ‘mwc’.
%
% The normalised, unmodulated images will be resampled to MNI space. These will 
% be written with the prefix ‘crwc’.
%
% =========================================================================

%-----------------------------------------------------------------------
% Job configuration created by cfg_util
%-----------------------------------------------------------------------
spm('defaults','fmri');
spm_jobman('initcfg');

matlabbatch{1}.spm.tools.preproc8.channel.vols = {[t1file,',1']};
matlabbatch{1}.spm.tools.preproc8.channel.biasreg = 0.0001;
matlabbatch{1}.spm.tools.preproc8.channel.biasfwhm = 60;
matlabbatch{1}.spm.tools.preproc8.channel.write = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).tpm = {[spmdir,'toolbox/Seg/TPM.nii,1']};
matlabbatch{1}.spm.tools.preproc8.tissue(1).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(1).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(1).warped = [1 1];
matlabbatch{1}.spm.tools.preproc8.tissue(2).tpm = {[spmdir,'toolbox/Seg/TPM.nii,2']};
matlabbatch{1}.spm.tools.preproc8.tissue(2).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(2).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(2).warped = [1 1];
matlabbatch{1}.spm.tools.preproc8.tissue(3).tpm = {[spmdir,'toolbox/Seg/TPM.nii,3']};
matlabbatch{1}.spm.tools.preproc8.tissue(3).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(3).native = [1 0];
matlabbatch{1}.spm.tools.preproc8.tissue(3).warped = [1 1];
matlabbatch{1}.spm.tools.preproc8.tissue(4).tpm = {[spmdir,'toolbox/Seg/TPM.nii,4']};
matlabbatch{1}.spm.tools.preproc8.tissue(4).ngaus = 3;
matlabbatch{1}.spm.tools.preproc8.tissue(4).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(4).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).tpm = {[spmdir,'toolbox/Seg/TPM.nii,5']};
matlabbatch{1}.spm.tools.preproc8.tissue(5).ngaus = 4;
matlabbatch{1}.spm.tools.preproc8.tissue(5).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(5).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).tpm = {[spmdir,'toolbox/Seg/TPM.nii,6']};
matlabbatch{1}.spm.tools.preproc8.tissue(6).ngaus = 2;
matlabbatch{1}.spm.tools.preproc8.tissue(6).native = [0 0];
matlabbatch{1}.spm.tools.preproc8.tissue(6).warped = [0 0];
matlabbatch{1}.spm.tools.preproc8.warp.mrf = 0;
matlabbatch{1}.spm.tools.preproc8.warp.reg = 4;
matlabbatch{1}.spm.tools.preproc8.warp.affreg = 'mni';
matlabbatch{1}.spm.tools.preproc8.warp.samp = 3;
matlabbatch{1}.spm.tools.preproc8.warp.write = [0 0];

% breakdown filename
[t1dir t1name t1suff] = fileparts(t1file);

% resample to MNI dimensions
matlabbatch{2}.spm.spatial.coreg.write.ref = {[spmdir,'templates/T1.nii,1']};
matlabbatch{2}.spm.spatial.coreg.write.source = {
                                                 [t1dir,'/wc1',t1name,t1suff,',1']
                                                 [t1dir,'/wc2',t1name,t1suff,',1']
                                                 [t1dir,'/wc3',t1name,t1suff,',1']
                                                 };
matlabbatch{2}.spm.spatial.coreg.write.roptions.interp = 1;
matlabbatch{2}.spm.spatial.coreg.write.roptions.wrap = [0 0 0];
matlabbatch{2}.spm.spatial.coreg.write.roptions.mask = 0;
matlabbatch{2}.spm.spatial.coreg.write.roptions.prefix = 'cr';

spm_jobman('run',matlabbatch);
