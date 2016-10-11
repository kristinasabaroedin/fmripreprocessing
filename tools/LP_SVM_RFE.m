	data = rand(60,100);
	labels = randi([0 1],60,1);
	labels(labels == 0) = -1;
function [r,r_srt] = LP_SVM_RFE(data,labels)
	% Input
	% data: subject by feature matrix
	% labels: vector of 1 and -1 denoting the allocation of subject's to one of two groups
	% Note, labels MUST be 1 and -1
	% 
	% Output
	% r: feature removal order
		% If numFeatures = 100
		% Then, r(1) = 61, means that feature there are 60 features more important than feature #1
	% 
	% r_srt: ranked feature list
		% the 1st element in r_srt is the most important single feature for SVM classifying.
		% r_srt(1) = 86, means that the most important feature is feature #86
		% r_srt(2) = 11, means that feature #11 is the second most important feature


	% Find number of features in data
	numFeatures = size(data,2);

	% Initialise feature ranking vector
	r = zeros(1,numFeatures);

	% Loop over features
	for i = 1:numFeatures

	    % Find features that are yet to be eliminated
	    s = find(r == 0);

	    % Check if there any features left to eliminate.
	    if length(s) > 1;

	    	% Pull out data from remaining features
	    	i_data = data(:,s);

	    	% Mdl = fitcsvm(i_data,labels);
	    	Mdl = svmtrain(i_data,labels);

	    	i_sv = Mdl.SupportVectors;
	    	i_alpha = Mdl.Alpha;
	    	i_labels = labels(Mdl.SupportVectorIndices);

	    	% FYI:
	    	i_weights = sum((i_alpha.*i_labels)'*i_sv);

	    	% Calculate weights for current features
	  %   	i_weights = zeros(1,numFeatures+1-i);
	  %   	for j = 1:size(i_sv,1); % loop over support vectors
			% 	i_weights = i_weights + i_alpha(j) * i_labels(j) * i_sv(j,:);
			% end

			% Find lowest absolute weight
	        [~, minweight] = min(abs(i_weights));
	        % Define weight to be removed
	        idx_remove = s(minweight);
	    else
	    	% If there is only one feature left, then use that as the idx
	    	idx_remove = s(1);
	    end

	    % Store the rank of the feature to be removed
	    r(1,idx_remove) = numFeatures+1-i;
	end

	% sort r into a ranked list of features
	[x, r_srt] = sort(r);
end