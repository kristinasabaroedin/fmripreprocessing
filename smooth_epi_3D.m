function [] = smooth_epi(epifile,kernel,N)

% smooth_epi(epifile,kernel,N)
%
% This function will perform spatial smoothing of an input fMRI timeseries.
%
% -------
% INPUTS:
% -------
%
% epifile   - a string containing the full path and filename of a 4D input 
%            EPI file. The images should have been co-registered to the T1. 
%            That is, they should be in T1 space; e.g.,
%            '/path/to/epi/4d.nii'
%
% kernel    - scalar value defining the size of the smoothing kernel in mm.
%
% N       - Lenth of time series (no of volumes in EPI 4D file).
%
% -------
% OUTPUTS:
% -------
%
% smoothed files will be output with the prefix 's'.
%
% =========================================================================

%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 3944 $)
%-----------------------------------------------------------------------

spm_jobman('initcfg')

for i = 1:N

    if i < 10
        ZeroPad = '0000';
    elseif i >= 10 & i < 100
        ZeroPad = '000';
    elseif i >= 100
        ZeroPad = '00';
    end
   
	matlabbatch{1}.spm.spatial.smooth.data{i,1} = [epifile,'_',ZeroPad,num2str(i),'.nii'];
end

matlabbatch{1}.spm.spatial.smooth.fwhm = repmat(kernel,1,3);  % define smoothing kernel based on input kernel value
matlabbatch{1}.spm.spatial.smooth.dtype = 0;
matlabbatch{1}.spm.spatial.smooth.im = 0;
matlabbatch{1}.spm.spatial.smooth.prefix = 's';

spm_jobman('run',matlabbatch);