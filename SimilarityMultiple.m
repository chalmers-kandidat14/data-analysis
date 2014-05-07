%%
%Imports chain data from csv files output
chainCsvFiles1 = dir('Chain2/chains/*.csv');
chainCsvFiles2 = dir('Chain3/chains/*.csv');
n = length(chainCsvFiles1);
chains = cell(n*2,1);

for k = 1:n
  filename = chainCsvFiles1(k).name;
  chains{k} = importdata(strcat('Chain2/chains/', filename),' ');
end

for k = 1:n
  filename = chainCsvFiles2(k).name;
  chains{k+n} = importdata(strcat('Chain3/chains/', filename),' ');
end
%%
%Imports graph data from csv files output, this is not needed as of now
graphCsvFiles1 = dir('Chain2/graphs/*.csv');
graphCsvFiles2 = dir('Chain3/graphs/*.csv');
graphs = cell(2*n,1);

for k = 1:n
  filename = graphCsvFiles1(k).name;
  graphs{k} = importdata(strcat('Chain2/graphs/', filename),' ');
end

for k = 1:n
  filename = graphCsvFiles2(k).name;
  graphs{k+n} = importdata(strcat('Chain3/graphs/', filename),' ');
end

%%
% Import energies
energy1 = importdata('Chain2/energies.csv');
energy2 = importdata('Chain3/energies.csv');

%%
%Calculates distance matrix by the root mean square measure of chain
%similarity, this takes some time and it is advisible to save the data to a
%dat file after this step

rmsSquareForm = zeros(2*n,2*n);
for i = 1:2*n
  for j = i+1:2*n
    [~,~, d] = Kabsch(chains{i},chains{j});
    rmsSquareForm(i,j) = d;
    rmsSquareForm(j,i) = d;
  end
  i
end
rmsVectorForm = squareform(rmsSquareForm);

%%
% Calculates the distance matrix by the graph method

graphSquareForm = zeros(2*n,2*n);
for i = 1:2*n
  for j = i+1:2*n
     dist = GraphDist(graphs{i},graphs{j});
     graphSquareForm(i,j) = dist;
     graphSquareForm(j,i) = dist;
  end
  i
end
graphVectorForm = squareform(graphSquareForm);

%%
%Creates the cluster tree and calculates the cophenet coefficient, which is
%a measure of how good the clustering is. 1 is best, 0 is worst.
clusterTreeMethod = 'average'; %'single', 'average' or 'complete'

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
numberOfClusters = 2;

clusters = cluster(clusterTree, 'maxclust', numberOfClusters);
chainColors = [ones(n,1) ; ones(n,1) + 1];
energy = -1*[energy1 ; energy2];

figure(3);
[Y, eigs] = cmdscale(graphSquareForm);

%colormap('copper')
scatter(Y(:,1),Y(:,2), 20, clusters)