function [] = first_level_contrasts(useTriStri,matfile)
%
%
% ------
% INPUTS
% ------
%
% matfile 	- path and file name of first level estimation output .mat file (typically SPM.mat)
%
% -------
% OUTPUTS
% -------
%
% SPM.mat file containing contrats for ROIs for 1st level analysis. located in spmdir.
%
% =========================================================================

spm('defaults','fmri');
spm_jobman('initcfg')
% run contrasts
if useTriStri == 1
	fprintf(1, 'Estimating first level contrasts using TriStri\n');
	matlabbatch{1}.spm.stats.con.spmmat = {matfile};
	matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Sensorimotor';
	matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 0 0 0];
	matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
	matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = 'Associative';
	matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [0 1 0 0];
	matlabbatch{1}.spm.stats.con.consess{2}.tcon.sessrep = 'none';
	matlabbatch{1}.spm.stats.con.consess{3}.tcon.name = 'Ventral';
	matlabbatch{1}.spm.stats.con.consess{3}.tcon.convec = [0 0 1 0];
	matlabbatch{1}.spm.stats.con.consess{3}.tcon.sessrep = 'none';
	matlabbatch{1}.spm.stats.con.delete = 0;
elseif useTriStri == 0
	fprintf(1, 'Estimating first level contrasts using sphere ROIs\n');
	matlabbatch{1}.spm.stats.con.spmmat = {matfile};
	matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = 'Ventral';
	matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [1 0];
	matlabbatch{1}.spm.stats.con.consess{1}.tcon.sessrep = 'none';
	matlabbatch{1}.spm.stats.con.delete = 0;
end

spm_jobman('run',matlabbatch);