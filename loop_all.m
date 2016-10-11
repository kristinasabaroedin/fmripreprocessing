clear all
sublist='~/Dropbox/scripts/projects/OCDPG/sublists/scratch.txt';
fileID = fopen(sublist);
subs = textscan(fileID,'%s');
subs = subs{1};

for i = 1:length(subs)
	fprintf(1,'Processing subject %s\n',subs{i})
	run_all_first(subs{i})
end
