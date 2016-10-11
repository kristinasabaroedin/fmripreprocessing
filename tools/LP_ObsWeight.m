%% ObsWeight: This function generates observation weights for use with SVM classification in cases where there is a class imbalance
function [Weights] = ObsWeight(Labels)

	ClassNames = unique(Labels);
	numClasses = size(ClassNames,1);

	Weights = zeros(size(Labels));

	for x = 1:numClasses
		y = ClassNames(x);
		Weights(Labels == y) = 1/(sum(Labels == y)/length(Labels));
	end
end