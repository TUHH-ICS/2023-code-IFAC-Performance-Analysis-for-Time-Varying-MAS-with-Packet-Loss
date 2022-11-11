function G = circular_graph(nvert, nedge, directed)
%CIRCULAR_GRAPH Generates a Matlab graph object that represents a circular
%graph with a given numer of vertices.
%   This function genenerates a Matlab graph object that contains a
%   circular graph with the given amount of vertices. It can be configured
%   if the connections should be directed or not and to how many of the
%   next vertices the connection should be established.
%
%   Arguments:
%       nvert    -> Number of vertices
%       nedge    -> Number of forward edges per vertex
%       directed -> [optional] Directed edges or not

if nargin <= 1 || isempty(nedge)
    nedge = 1;
elseif nedge >= nvert
    error('So many neighbours are not possible!')
end
if nargin <= 2
    directed = true;
end

% Define closure that keeps the index in the ring [1, nvert]
ring = @(i) mod(i-1, nvert) + 1;

%% Generate graph structure
A = zeros(nvert);
for i = 1:nvert
    A(i, ring(i+(1:nedge))) = 1;
    
    % If not direction, also add the other direction
    if ~directed
        A(i, ring(i-(1:nedge))) = 1;
    end
end

% G needs only to be a digraph if it is directed
if directed
    G = digraph(A);
else
    G = graph(A);
end
end
