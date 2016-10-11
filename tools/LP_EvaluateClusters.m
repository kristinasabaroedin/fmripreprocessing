function LP_EvaluateClusters(dataMatrix)

krange = 1:10;
numReplicates = 10;
dmeth = 'sqEuclidean';

fprintf(1,'Evaluating optimal clustering solution...\n');

%-------------------------------------------------------------------------------
% Do the (kmeans) clustering across k range
rng(0);
clusterings = zeros(size(dataMatrix,1),length(krange));
mykmeans = @(X,K) kmeans(X,K,'distance',dmeth,...
                            'replicates',numReplicates, ...
                            'start','sample');

for i = 1:length(krange)
    clusterings(:,i) = mykmeans(dataMatrix,krange(i));
end

% ------------------------------------------------------------------------------
% Optimal k (number of clusters) according to different criteria:
% ------------------------------------------------------------------------------
eva_CH = evalclusters(dataMatrix,clusterings,'CalinskiHarabasz');
eva_DB = evalclusters(dataMatrix,clusterings,'DaviesBouldin');
eva_Gap = evalclusters(dataMatrix,mykmeans,'gap','B',10,'KList',krange,...
                'ReferenceDistribution','uniform',...
                'SearchMethod','globalMaxSE');
eva_Sil = evalclusters(dataMatrix,clusterings,'silhouette');


% ------------------------------------------------------------------------------
% Visualize results:
% ------------------------------------------------------------------------------
f = figure('color','w');
krange = 1:10;
disp(eva_Gap)
subplot(2,3,1)
plot(eva_Gap.LogW,'o-k'); xlabel('k'); ylabel('Observed and Expected log','interpreter','none')
hold on
plot(eva_Gap.ExpectedLogW,'x-k')
subplot(2,3,2)
plot(eva_Gap)
title(sprintf('Gap statistic, optimal k = %u',eva_Gap.OptimalK))
subplot(2,3,3)
plot(eva_Sil)
title(sprintf('Silhouette statistic, optimal k = %u',eva_Sil.OptimalK))
subplot(2,3,4)
plot(eva_CH)
title(sprintf('Calinski-Harabasz criterion, optimal k = %u',eva_CH.OptimalK))
subplot(2,3,5)
plot(eva_DB)
title(sprintf('Davies Bouldin criterion, optimal k = %u',eva_DB.OptimalK))

end
