%%
%Chains
%chainfolders = {'Chain1', 'Chain2', 'Chain3','Chain4', 'Chain5', 'Chain6','Chain7', 'Chain8', 'Chain9'}
chainfolders = {'Chain2Cubic'}
nochains = length(chainfolders);
importsperchain = 1000;

%%
%Imports chain data from csv files output

n = nochains*importsperchain;S
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
      a = importdata(filepath, ' ');
      graphs{index}=a+triu(a')+triu(a)';
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
     dist = GraphDist(graphs2{i},graphs2{j});
     squareForm(i,j) = dist;
     squareForm(j,i) = dist;
  end
end
 
vectorForm = squareform(squareForm);

%%
%Calculates the medioid
meandists = mean(squareForm, 1);
[~, medioid] = min(meandists);
medioid

%% Imports the energies
energies = importdata([chainfolders{1} '/energies.csv']);

%% Calculate what neighbours are common
commonneigh = graphs{1};
for i = 2:n
    commonneigh = commonneigh + graphs{i};
end
commonneigh = commonneigh / n;
sortedcommonneigh = sort(reshape(commonneigh, [],1));

%% Calculate what neighbours are common in hpstruct
hdir=dir('graph-*.csv');
hgraphs=cell(13,1);
for i=1:13
  filename=hdir(i).name;
  a=importdata(filename, ' ');
  hgraphs{i}=a+triu(a')+triu(a)';
end

hcommonneigh = hgraphs{1};
for i = 2:13
    hcommonneigh = hcommonneigh + hgraphs{i};
end
hcommonneigh = hcommonneigh / 13;
hsortedcommonneigh = sort(reshape(hcommonneigh, [],1));

%% Calculates total amount of neighbours
noneigh = ones(1, n);
for i = 1:n
    noneigh(i)=sum(sum(graphs{i}));
end
avrgneighbs = mean(noneigh)

%%
%Creates the cluster tree and calculates the cophenet coefficient, which is
%a measure of how good the clustering is. 1 is best, 0 is worst.
clusterTreeMethod = 'average'; %'single', 'average' or 'complete'

clusterTree = linkage(vectorForm, clusterTreeMethod);
'The cophenet coefficien of the cluster tree is:'
cophenet(clusterTree, vectorForm)

%%
%Creates a dendrogram of the cluster tree

figure(1)
[H, T, outperm] = dendrogram(clusterTree);

%%
%Calculates and displays a multi-dimentional scaling visualization of the
%data with clusters colorized.
numberOfClusters = 5;

clusters = cluster(clusterTree, 'maxclust', numberOfClusters);
figure(2);
chainColor = [ones(1000,1) ; ones(13,1)+1];

[Y, eigs] = cmdscale(squareForm);
scatter(Y(:,1),Y(:,2), 50, chainColor)

