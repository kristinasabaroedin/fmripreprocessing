function [fd,mask,delta_mov] = scrub_fd(mov, thr, head)

	% This function differentiates head movement parameters across frames and 
	% computes a measure of total framewise displacement over time, as detailed 
	% Power et al. (2012) NeuroImage. It also identifies time points exceeding 
	% a specific displacement thresholds which may be subsequently excluded
	% ('scrubbed') from the time series.
	%
	% ------
	% INPUTS
	% ------
	% mov       - an N x 6 matrix containing 6 movement parameters and where N
	%           = length of the time series.
	% thr       - threshold for identifying frames that should be excluded
	%           (default = .5, as in Power et al.)
	% head      - head radius (in mm) to use when converting radians to mm. default =
	%           50mm, as in Power et al.
	%
	% -------
	% OUTPUTS
	% -------
	% fd        - an N-1 length vector representing the total framewise
	%           displacement
	% mask      - indices of frames that have fd > thr
	% delta_mov - N-1 x 6 matrix of differentiated movement parameters
	%
	% Alex Fornito, Melbourne Neuropsychiatry Centre, 2012
	%__________________________________________________________________________

	% convert degrees to radians
	mov(:,4:6) = degtorad(mov(:,4:6));

	% convert radians to mm
	mov(:,4:6) = mov(:,4:6).*head;

	% differentiate movement parameters
	for i = 1:size(mov,2)
	    delta_mov(:,i) = diff(mov(:,i));
	end

	% compute total framewise displacement
	fd = sum(abs(delta_mov'));
	   
	% find frames exceeding threshold
	mask = find(fd>thr);
	mask = mask+1;  % add one because second value in pair should be removed from original time series.

end      