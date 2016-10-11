%% ConfAcc: 
function [Accuracy] = ConfAcc(ConfMat)
	for i = 1:size(ConfMat,1)
		Accuracy(i) = ConfMat(i,i)/sum(ConfMat(:,i));
	end
end
