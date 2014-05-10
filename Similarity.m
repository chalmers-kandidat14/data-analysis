%%
%Chains
%chainfolders = {'Chain1', 'Chain2', 'Chain3','Chain4', 'Chain5', 'Chain6','Chain7', 'Chain8', 'Chain9'}
chainfolders = {'Chain8Cubic','Chain9Cubic'}
nochains = length(chainfolders);
importsperchain = 1000;

%%
%Imports chain data from csv files output

n = nochains*importsperchain;
chains = cell(n,1);
for i = 1:nochains;
    chainCsvFiles = dir([chainfolders{i}, '/chains/*.csv'])
    for k = 1:importsperchain
      filename = chainCsvFiles(k).name;
      index = (i-1)*importsperchain+k;
      filepath = [chainfolders{i}, '/chains/', filename];
      chains{index} = importdata(filepath, ' ');
    end
end
%%
%Imports graph data from csv files output, this is not needed as of now
n = nochains*importsperchain;
graphs = cell(n,1);
for i = 1:nochains;
    i
    graphCsvFiles = dir([chainfolders{i}, '/graphs/*.csv']);
    for k = 1:importsperchain
      filename = graphCsvFiles(k).name;
      index = (i-1)*importsperchain+k;
      filepath = [chainfolders{i}, '/graphs/', filename];
      graphs{index} = importdata(filepath, ' ');
    end
end
%%
%Calculates distance matrix by the root mean square measure of chain
%similarity, this takes some time and it is advisible to save the data to a
%dat file after this step

squareForm = zeros(n,n);
for i = 1:n
  for j = i+1:n
    [~,~, d] = Kabsch(chains{i},chains{j});
    squareForm(i,j) = d;
    squareForm(j,i) = d;
  end
  i
end
vectorForm = squareform(squareForm);

%%
% Calculates the distance matrix by the graph method

squareForm = zeros(n,n);
for i = 1:n
  for j = i+1:n
     dist = GraphDist(graphs{i},graphs{j});
     squareForm(i,j) = dist;
     squareForm(j,i) = dist;
  end
end
vectorForm = squareform(squareForm);

%%
%Creates the cluster tree and calculates the cophenet coefficient, which is
%a measure of how good the clustering is. 1 is best, 0 is worst.
clusterTreeMethod = 'average'; %'simple', 'average' or 'complete'

clusterTree = linkage(graphVectorForm, clusterTreeMethod);
'The cophenet coefficien of the cluster tree is:'
cophenet(clusterTree, graphVectorForm)

%%
%Creates a dendrogram of the cluster tree

figure(1)
[H, T, outperm] = dendrogram(clusterTree);

%%
%Calculates and displays a multi-dimentional scaling visualization of the
%data with clusters colorized.
numberOfClusters = 9;

clusters = cluster(clusterTree, 'maxclust', numberOfClusters);
figure(3);

[Y, eigs] = cmdscale(graphSquareForm);
scatter(Y(:,1),Y(:,2), 50, clusters)

