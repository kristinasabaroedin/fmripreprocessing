function [] = first_level(resdir,units,TR,data,N,matfile,mask)

% first_level(spmdir,timing,TR,data,N,matfile)
%
% This script will run 1st level specification and estimation for a
% seed-based analysis. That is, only 'multiple regressors' are
% specified--nothing that requires haemodynamic convolution.
%
% ------
% INPUTS
% ------
%
% resdir    - string representing path where results will be written to.
%
% timing    - timing unit for analysis. Should be either 'scans' or
%           'seconds'.
%
% TR        - repetition time (TR) of data
%
% data      - path and filename of epi 4d file on which the analysis will
%           be run (e.g., smoothed and normalized epis)
%
% N         - Length of time series (no. of vols)
%
% matfile   - path and filename of .mat file containing regressors for the analysis. Should
%           contain and N*K matrix, where K is the number of regressors.
%           The matrix must be called 'R' as per SPM requirements.
%
% mask      - optional string indicating path and filename of an explicit
%           mask image (e.g., skull-stripped t1) to constrain analysis. If
%           no mask, set to ''.
%
% -------
% OUTPUTS
% -------
%
% SPM.mat file contating results of 1st level analysis. located in spmdir.
%
% =========================================================================

% resdir = [rawdir,'spm/'];
% units = 'scans';
% TR = 2;
% data = [rawdir,'squared.nii'];
% N = 234;
% matfile = [rawdir,'spm_regs'];



%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 3944 $)
%-----------------------------------------------------------------------

spm('defaults','fmri');
spm_jobman('initcfg')

for i = 1:N

    if i < 10
        ZeroPad = '0000';
    elseif i >= 10 & i < 100
        ZeroPad = '000';
    elseif i >= 100
        ZeroPad = '00';
    end
    
    matlabbatch{1}.spm.stats.fmri_spec.sess.scans{i,1} = [data(1:end-4),'_',ZeroPad,num2str(i),'.nii'];
    
end

matlabbatch{1}.spm.stats.fmri_spec.dir = {resdir};
matlabbatch{1}.spm.stats.fmri_spec.timing.units = units;
matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = 16;
matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = 1;
matlabbatch{1}.spm.stats.fmri_spec.sess.cond = struct('name', {}, 'onset', {}, 'duration', {}, 'tmod', {}, 'pmod', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi = {''};
matlabbatch{1}.spm.stats.fmri_spec.sess.regress = struct('name', {}, 'val', {});
matlabbatch{1}.spm.stats.fmri_spec.sess.multi_reg = {matfile};
matlabbatch{1}.spm.stats.fmri_spec.sess.hpf = -Inf;
matlabbatch{1}.spm.stats.fmri_spec.fact = struct('name', {}, 'levels', {});
matlabbatch{1}.spm.stats.fmri_spec.bases.hrf.derivs = [0 0];
matlabbatch{1}.spm.stats.fmri_spec.volt = 1;
matlabbatch{1}.spm.stats.fmri_spec.global = 'None';
matlabbatch{1}.spm.stats.fmri_spec.mask = {mask};
matlabbatch{1}.spm.stats.fmri_spec.cvi = 'none';%'AR(1)';

% run estimation
matlabbatch{2}.spm.stats.fmri_est.spmmat = {[resdir,'SPM.mat']};
matlabbatch{2}.spm.stats.fmri_est.method.Classical = 1;

spm_jobman('run',matlabbatch);