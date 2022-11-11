function L = laplace_matrix(G)
%LAPLACE_MATRIX Calculate the Laplacian of the given graph
%   This function is a replacement for the Matlab built-in laplacian()
%   function that also works for digraphs, unlike the Matlab version.

A = adjacency(G);
D = diag(sum(A,2));
L = D - A;
end
