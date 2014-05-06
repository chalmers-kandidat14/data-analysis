function [ d ] = GraphDist( A, B )
%GRAPHSIM Calculates the distance between two graphs
%   Detailed explanation goes here
  
d = sum(sum(A) + sum(B) - 2*sum(A.*B));

end