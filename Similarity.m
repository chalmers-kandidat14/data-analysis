%%
%Imports chain data from csv files output
chainCsvFiles = dir('chains/*.csv');
n = length(chainCsvFiles);
chains = cell(n,1);

for k = 1:n
  filename = chainCsvFiles(k).name;
  chains{k} = importdata(strcat('chains/', filename),' ');
end
%%
%Imports graph data from csv files output, this is not needed as of now
graphCsvFiles = dir('graphs/*.csv');
n = length(graphCsvFiles);
graphs = cell(n,1);

for k = 1:n
  filename = graphCsvFiles(k).name;
  graphs{k} = importdata(strcat('graphs/', filename),' ');
end
%%
%Calculates distance matrix by the root mean square measure of chain
%similarity, this takes some time and it is advisible to save the data to a
%dat file after this step

rmsSquareForm = zeros(n,n);
for i = 1:n
  for j = i+1:n
    [~,~, d] = Kabsch(chains{i},chains{j});
    rmsSquareForm(i,j) = d;
    rmsSquareForm(j,i) = d;
  end
  i
end
rmsVectorForm = squareform(rmsSquareForm);
%%
%Creates the cluster tree and calculates the cophenet coefficient, which is
%a measure of how good the clustering is. 1 is best, 0 is worst.
clusterTreeMethod = 'average'; %'simple', 'average' or 'complete'

clusterTree = linkage(rmsVectorForm, clusterTreeMethod);
'The cophenet coefficien of the cluster tree is:'
cophenet(clusterTree, rmsVectorForm)
%%
%Creates a dendrogram of the cluster tree

figure(1)
[H, T, outperm] = dendrogram(clusterTree);

%%
%Calculates and displays a multi-dimentional scaling visualization of the
%data with clusters colorized.
numberOfClusters = 10;

clusters = cluster(clusterTree, 'maxclust', numberOfClusters);
figure(3);
[Y, eigs] = cmdscale(rmsSquareForm);
scatter(Y(:,1),Y(:,2), 50, clusters)