% SVM_ClassAcc: This computes robust estimates of classification accuracy across a range of PCs
% using N repeats of k-fold cross validation
function [ClassAcc_Iter,ClassAcc,ClassAccStd,ConfMat_Iter,ClustAcc_Iter,ClustAcc,ClustAccStd] = SVM_ClassAcc(DataIn,Labels,iterations,t,switchFeatures,Weights)
	% Inputs
	%
	% DataIn 				- data for classification. Rows correspond to observations (e.g., participants) and columns to features
	% Labels 				- a vector of the groupings of the DataIn rows. This serves as the ground truth for classification.
	% iterations 			- the number of iterations for SVM repeats. More repeats makes for more robust estimates of classification metrics.
	% t 					- SVM template
	% 							e.g., 	t = templateSVM('Standardize',1,'KernelFunction','linear');
	% switchFeatures		- Switch that determines how the features are input to the SVM.
	% 							1 = Use all features.
	% 							2 = use single features, entered one at a time.
	% 							3 = input features as incremental groups.
	% Weights (optional)	- observation weights for SVM
	% ------------------------------------------------------------------------------
	numSamples = size(DataIn,1);
	Features = 1:size(DataIn,2);

	if nargin < 6
		Weights = ones(1,numSamples);
	end

	% ------------------------------------------------------------------------------
	% Use ALL features concurrently
	% ------------------------------------------------------------------------------
	if switchFeatures == 1
		% Storage variables
		ClassAcc_Iter = zeros(1,iterations);
		ConfMat_Iter = cell(1,iterations);
		ClustAcc_Iter = cell(1,iterations);

		rng('default');
		fprintf(1,'Performing SVM w/ all features\n')
		numFeatures = Features(end);

		for i = 1:iterations
			fprintf(1,'...%u \n',i)
	
			% Fit the SVM model
			Mdl = fitcecoc(DataIn,Labels,'Learners',t,'Coding','onevsall','Weights',Weights);
				
			% Cross validate
			CVMdl = crossval(Mdl);
			oofLabel(:,i) = kfoldPredict(CVMdl);

			% compute classification accuracy for each iteration
			ClassAcc_Iter(i) = sum(oofLabel(:,i) == Labels)/length(Labels);
			
			% compute confusion on a per cluster basis - allows for per cluster accuracy
			ConfMat_Iter{i} = confusionmat(oofLabel(:,i),Labels);
			ClustAcc_Iter{i} = LP_ConfAcc(ConfMat_Iter{i});
		end
		ClustAcc = []; ClustAccStd = [];
		% compute SVM per clust overall accuracy
		ClustAcc = mean(cat(1,ClustAcc_Iter{1,:}));
		ClustAccStd = std(cat(1,ClustAcc_Iter{1,:}));
	% ------------------------------------------------------------------------------
	% Use either individual features or binned features
	% ------------------------------------------------------------------------------
	elseif switchFeatures ~= 1
		% Storage variables
		ClassAcc_Iter = zeros(length(Features),iterations);
		ConfMat_Iter = cell(length(Features),iterations);
		ClustAcc_Iter = cell(length(Features),iterations);

		for j = 1:length(Features) % loop over number of PCs
			rng('default');
			numFeatures = Features(j);
			if switchFeatures == 2
				fprintf(1,'Performing SVM w/ Feature %u\n',j)
			elseif switchFeatures == 3
				fprintf(1,'Performing SVM w/ %u Features\n',numFeatures)
			end			

			oofLabel = zeros(numSamples,iterations);
			
			for i = 1:iterations
				fprintf(1,'...%u \n',i)
			
				% Fit the SVM model
				if switchFeatures == 2
					Mdl = fitcecoc(DataIn(:,j),Labels,'Learners',t,'Coding','onevsall','Weights',Weights);
				elseif switchFeatures == 3
					Mdl = fitcecoc(DataIn(:,1:numFeatures),Labels,'Learners',t,'Coding','onevsall','Weights',Weights);
				end

				% Cross validate
				CVMdl = crossval(Mdl);
				oofLabel(:,i) = kfoldPredict(CVMdl);

				% compute classification accuracy for each iteration
				ClassAcc_Iter(j,i) = sum(oofLabel(:,i) == Labels)/length(Labels);
				
				% compute confusion on a per cluster basis - allows for per cluster accuracy
				ConfMat_Iter{j,i} = confusionmat(oofLabel(:,i),Labels);
				ClustAcc_Iter{j,i} = LP_ConfAcc(ConfMat_Iter{j,i});
			end
		end
	
		ClustAcc = []; ClustAccStd = [];
		% compute SVM per clust overall accuracy
		for i = 1:length(Features)
			ClustAcc(i,:) = mean(cat(1,ClustAcc_Iter{i,:}));
			ClustAccStd(i,:) = std(cat(1,ClustAcc_Iter{i,:}));
		end
	end
	% compute SVM overall accuracy
	ClassAcc = mean(ClassAcc_Iter,2);
	ClassAccStd = std(ClassAcc_Iter,0,2);
end