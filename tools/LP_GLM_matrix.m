%% GLM_matrix: This function takes a cell with subject IDs and prepares a groupwise GLM matrix for use with FSL's ranomise
function [matrix] = GLM_matrix(subs)
	matrix = [];
	for i = 1:length(subs)
		x = ones(length(subs{i}),1);
		x = padarray(x,[0 length(subs)-i],0,'pos');
		x = padarray(x,[0 i-1],0,'pre');

		matrix = [matrix;x];
	end
end