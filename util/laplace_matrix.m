%---------------------------------------------------------------------------------------------------
% For Paper
% "Robust Performance Analysis for Time-Varying Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function L = laplace_matrix(G)
%LAPLACE_MATRIX Calculate the Laplacian of the given graph
%   This function is a replacement for the Matlab built-in laplacian()
%   function that also works for digraphs, unlike the Matlab version.

A = adjacency(G);
D = diag(sum(A,2));
L = D - A;
end
