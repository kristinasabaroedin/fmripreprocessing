function [] = prepro_reports(mni, t1, t12mni, wm2mni, csf2mni, epi2mni) 

% [] = prepro_reports(mni, t1, t1_brain, crt1, t12mni, epi2t1, epi2mni, epi_brain) 
%
% This function generates a series of reports displaying pre-processing
% results for quality control.
%
% ------
% INPUTS
% ------
% NB: All inputs are strings
%
% mni        - path and filename of template used to register to MNI space;
%              e.g., '/usr/loca/spm8/templates/T1.nii'
%
% t1         - path and filename of t1
%
% t1_brain   - path and filename of skull-stripped T1
%
% crt1       - path and filename of t1 co-registered to MNI
%
% t12mni     - path and filename of t1 normalized to MNI template
%
% epi2mni    - path and filename of epi spatially normalized to MNI
%              template.
%
% -------
% OUTPUTS
% -------
%
% A .ps file containing the following results:
%
% Check reg:    norm t1 vs MNI
%               norm epi vs MNI
%               norm wm mask vs norm epi
%               norm csf mask vs norm epi
%
%
% =========================================================================


%-----------------------------------------------------------------------
% Job configuration created by cfg_util (rev $Rev: 3944 $)
%-----------------------------------------------------------------------
spm('defaults','fmri');
spm_jobman('initcfg')

matlabbatch{1}.spm.util.checkreg.data = {
                                         [t12mni,',1']
                                         [mni,',1']
                                         };

matlabbatch{2}.spm.util.print.fname = 't1-mni_v_mni';
matlabbatch{2}.spm.util.print.fig.figname = 'Graphics';
matlabbatch{2}.spm.util.print.opts.opt = {
                                          '-dpsc2'
                                          '-append'
                                          }';
matlabbatch{2}.spm.util.print.opts.append = true;
matlabbatch{2}.spm.util.print.opts.ext = '.ps';

matlabbatch{3}.spm.util.checkreg.data = {
                                         [epi2mni,',1']
                                         [mni,',1']
                                         };

matlabbatch{4}.spm.util.print.fname = 'epi-mni_v_mni';
matlabbatch{4}.spm.util.print.fig.figname = 'Graphics';
matlabbatch{4}.spm.util.print.opts.opt = {
                                          '-dpsc2'
                                          '-append'
                                          }';
matlabbatch{4}.spm.util.print.opts.append = true;
matlabbatch{4}.spm.util.print.opts.ext = '.ps';

matlabbatch{5}.spm.util.checkreg.data = {
                                         [wm2mni,',1']
                                         [epi2mni,',1']
                                         };
                                     
matlabbatch{6}.spm.util.print.fname = 'wm-mni_v_epi-mni';
matlabbatch{6}.spm.util.print.fig.figname = 'Graphics';
matlabbatch{6}.spm.util.print.opts.opt = {
                                          '-dpsc2'
                                          '-append'
                                          }';
matlabbatch{6}.spm.util.print.opts.append = true;
matlabbatch{6}.spm.util.print.opts.ext = '.ps';   

matlabbatch{7}.spm.util.checkreg.data = {
                                         [csf2mni,',1']
                                         [epi2mni,',1']
                                         };
                                     
matlabbatch{8}.spm.util.print.fname = 'csf-mni_v_epi-mni';
matlabbatch{8}.spm.util.print.fig.figname = 'Graphics';
matlabbatch{8}.spm.util.print.opts.opt = {
                                          '-dpsc2'
                                          '-append'
                                          }';
matlabbatch{8}.spm.util.print.opts.append = true;
matlabbatch{8}.spm.util.print.opts.ext = '.ps'; 

spm_jobman('run',matlabbatch)
