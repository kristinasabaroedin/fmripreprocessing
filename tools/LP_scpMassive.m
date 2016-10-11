clear all

Project = 'NAMIC';
INdir = ['lindenmp@m2.massive.org.au:/gpfs/M2Home/lindenmp/Monash076_scratch/rannee/',Project,'/'];
OUTdir = ['/media/lindenmp/SG8_4TB_2/ResProjects/',Project,'/'];

sublist = [OUTdir,'SubjectList.txt'];
fileID = fopen(sublist);
subs = textscan(fileID,'%s');
subs = subs{1};

for i = 1:length(subs)
	fprintf(1,'Copying subject %s\n',subs{i})
	
	% ------------------------------------------------------------------------------
	% T1 data
	% ------------------------------------------------------------------------------
	fprintf(1, 'Copying T1 data...');

	t1dir = [OUTdir,subs{i},'/anat_1/'];
	if exist(t1dir) == 0
		mkdir(t1dir)
	elseif exist(t1dir) == 7
		rmdir(t1dir,'s')
		mkdir(t1dir)
	end

	system(['scp ',INdir,'/',subs{i},'/anat_1/mprage.nii ',t1dir])

	fprintf(1, 'done\n');

	% ------------------------------------------------------------------------------
	% rfMRI data
	% ------------------------------------------------------------------------------
	fprintf(1, 'Copying rfMRI data...');
	
	rfMRIdir = [OUTdir,subs{i},'/rest_1/'];
	if exist(rfMRIdir) == 0
		mkdir(rfMRIdir)
	elseif exist(rfMRIdir) == 7
		rmdir(rfMRIdir,'s')
		mkdir(rfMRIdir)
	end

	system(['scp ',INdir,'/',subs{i},'/rest_1/rest.nii ',rfMRIdir])
	
	fprintf(1, 'done\n');

end
