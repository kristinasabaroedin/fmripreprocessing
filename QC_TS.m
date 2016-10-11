clear all
% function [] = QC_TS(subject,mov,rsdata,t1)

datadir = '/gpfs/M2Home/kristina_s/Monash076/Kristina/GenofCog/data/';
subject = '1008.2.57.034';

rsdata = 'srest_prepro.nii';
t1 = 't1_segmask.nii';

	% load t1 binary mask
	cd([datadir,subject,'/t1/'])
	[hdr_t1,t1] = read(t1);


	% load in resting state data
	cd([datadir,subject,'/rfMRI/'])
	% cd([datadir,subject,'/rfMRI/'])
	[hdr_ts,ts] = read(rsdata);

	% Reshape to 2d matrix, with voxels on one axis and time points on the other
	dim = size(ts);
	ts_2d = reshape(ts,dim(1)*dim(2)*dim(3),dim(4));
	dim = size(t1);
	t1_2d = reshape(t1,dim(1)*dim(2)*dim(3),1);

	% remove non brain voxels from ts
	ts_2d(t1_2d == 0,:) = [];
	
	% read in motion
    mfile = dir('rp*txt');
    mov = dlmread([rawdir,mfile(1).name]);
    mov = detrend(mov,'linear'); % detrend motion regressors

	% plot
	figure
	subplot(3,1,1)
	plot(mov(:,1));
	hold on
	plot(mov(:,2),'g');
	plot(mov(:,3),'r');
	title('translation','fontsize',15,'fontweight','bold')
	legend({'x','y','z'})
	ylabel('mm')
	xlabel('time (volumes)')
	xlim([1 size(ts,2)])

	subplot(3,1,2)
	plot(mov(:,4));
	hold on
	plot(mov(:,5),'g');
	plot(mov(:,6),'r');
	title('rotation','fontsize',15,'fontweight','bold');
	legend({'pitch','roll','yaw'})
	ylabel('degrees')
	xlabel('time (volumes)')
	xlim([1 size(ts,2)])

	subplot(3,1,3)
	imagesc(ts_2d)
	colormap([flipud(BF_getcmap('blues',10));1,1,1;BF_getcmap('reds',10)])
	% colorbar	
	temp = max(abs(caxis));
	caxis([-temp/3 temp/3])

	title('Voxel time series (RED = peaks, BLUE = troughs)','fontsize',15,'fontweight','bold');
	ylabel('Voxels')
	xlabel('time (volumes)')

	set(gcf,'color','white');

% end