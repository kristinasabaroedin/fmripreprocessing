clear all; close all
rng('default')
Plot = 1;

sublist='~/Dropbox/scripts/projects/OCDPG/sublists/SubjectIDs.txt';
fileID = fopen(sublist);
subs = textscan(fileID,'%s');
subs = subs{1};

% projdir = '/media/lindenmp/SG8_4TB/Research_Projects/OCDPG';
projdir = '/Volumes/SG8_4TB/Research_Projects/OCDPG';

workdir = [projdir,'/data'];

numParams = 6;

maxMov = [];
fd = cell(length(subs),1);
mask_count = zeros(length(subs),1);
mask = cell(length(subs),1);
delta_mov = cell(length(subs),1);

for i = 1:length(subs)
	fprintf(1, 'Processing %s\n', char(subs(i)));
	str = char(strcat(workdir,'/',subs(i),'/rfMRI'));
	cd(str)
	% ------------------------------------------------------------------------------
	% Load in noise_signals.txt
	% ------------------------------------------------------------------------------
    mov = dlmread('noise_signals.txt');

    % ------------------------------------------------------------------------------
    % Find max values for translation and rotation
    % ------------------------------------------------------------------------------
	maxMov(i,:) = max(abs(mov(:,1:numParams)));

	% ------------------------------------------------------------------------------
	% Compute Power (2012) scrub metric
	% ------------------------------------------------------------------------------
	[fd{i},mask{i},delta_mov{i}] = scrub_fd(mov(:,1:numParams),0.5,50);

	% ------------------------------------------------------------------------------
	% vector that denotes how many fd problems a participant has
	% ------------------------------------------------------------------------------
	if ~isempty(mask{i})
		mask_count(i) = numel(mask{i});
	else
		mask_count(i) = 0;
	end
		
end


% ------------------------------------------------------------------------------
% Save data
% ------------------------------------------------------------------------------
outdir='/Users/lindenmp/Dropbox/University/Monash/PhD/Research_Projects/OCD-PG';
cd(outdir)
save('fd.mat','subs','maxMov','delta_mov','fd','mask_count','mask')

% ------------------------------------------------------------------------------
% Plot
% ------------------------------------------------------------------------------
    if Plot
	    outdir = [projdir,'/rfMRI_fd'];
	    % ------------------------------------------------------------------------------
	    % Initialise output directory. If already exists, delete and re-initialise
	    % ------------------------------------------------------------------------------
	    if exist(outdir) == 0
	        fprintf(1,'Initialising outdir\n')
	        mkdir(outdir)
	    elseif exist(outdir) == 7
	        fprintf(1,'Cleaning and re-initialising outdir\n')
	        rmdir(outdir,'s')
	        mkdir(outdir)
	    end
	    cd(outdir)
		
		for i = 1:length(subs)
			if ~isempty(mask{i})
				close all
				plot(fd{i})
				title(subs(i))
				saveas(gcf,[char(subs(i)),'.pdf'],'pdf')
			end
		end
	end

