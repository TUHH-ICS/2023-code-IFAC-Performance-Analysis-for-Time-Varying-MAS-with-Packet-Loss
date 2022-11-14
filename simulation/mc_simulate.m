%---------------------------------------------------------------------------------------------------
% For Paper
% "Robust Performance Analysis for Time-Varying Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

function [h2, converged] = mc_simulate(graphs, pattern, simconf, netconf)
%MC_SIMULATE Perform a Monte-Carlo simulation of the H2-norm of the
%network of connected masses
%   This function repeatedly runs simulations of networks of masses trying
%   to maintain a formation under disturbances. By varying the properties
%   of the network, this can be exploited to test controllers under varying
%   circumstances.
%
%   Arguments:
%       graphs  -> Nominal communication graphs
%       pattern -> Switching pattern between the graphs
%       simconf -> Configuration for the simulation
%       netconf -> Might be struct array, configuration for network. Either
%                  SinrConfiguration or plain struct for Bernoulli
%   Returns:
%       h2        -> Mean H2-norm of the system for each type of network
%       converged -> True, if all simulations converged

% Check consistency between pattern and provided graphs
if min(pattern, [], 'all') < 1 || max(pattern, [], 'all') > length(graphs)
    error('The switching pattern is not consistent with the number of provided graphs')
end

% Extract problem data
[dim, N] = size(simconf.positions);
K = length(graphs);

% Calculate all required graph matrices
L0 = zeros(N,N,K);
A0 = zeros(N,N,K);
Ld = zeros(N*dim,N*dim,K);
for i = 1:length(graphs)
    if height(graphs{i}.Nodes) ~= N
        error('Graph %d does not fit the problem size', i)
    end

    L0(:,:,i) = full(laplace_matrix(graphs{i}));
    A0(:,:,i) = full(adjacency(graphs{i}));
    Ld(:,:,i) = kron(L0(:,:,i), eye(dim));
end
Rhat = kron(eye(N), sqrtm(simconf.R));
Pi   = kron(eye(N) - ones(N)/N, eye(dim));

% Allocate storage for the simulation parameters and results
n = 1:simconf.samples;
m = 1:N*dim;
[Net, Pat, M, ~] = ndgrid(netconf, 1:size(pattern,1), m, n);
perf = zeros(size(M));
conv = zeros(size(M), 'logical'); % Convergence flag

% Create progress meter
meter = ParforProgressMeter(numel(M), 0.002);
meter.start()

parfor i = 1:numel(M)
    % Initialize the network
    if isa(Net, 'SinrConfiguration')
        Network = SinrNetwork(Net(i));
    else
        Network = BernoulliNetwork(N, simconf.dT, dim,...
                                    Net(i).range, Net(i).p, Net(i).sym);
    end

    % Extract current pattern
    cur_pat = pattern(Pat(i),:);

    % Initialize the agents
    Agents = cell(N, 1);
    for j = 1:length(Agents)
        pos  = simconf.positions(:,j);
        dist = (mod(M(i)-1, dim)+1) * (mod(M(i)-1, N) == j-1);
        adj  = reshape(A0(j,:,:), N, K);
        Agents{j} = DisturbedAgent(Network.getId(), pos, pos, adj, cur_pat, dist);
    end
    Agents = [Agents{:}];

    % Initialize simulation
    sim  = SimulationManager(Network, Agents);
    idx  = Agents(1).pattern_idx;
    init = true;
    
    % Simulate!
    t = sim.step();
    while formation_error(Agents) > simconf.tol && t < simconf.Tf       
        % Evaluate formation error
        pos = [Agents.position] - [Agents.ref];
        err = Ld(:,:,cur_pat(idx))*pos(:);
        
        % Skip control input in very first iteration, because it only
        % contains the impulse input
        if init
            ctr  = 0;
            init = false;
        else
            u   = [Agents.u];
            ctr = Pi*Rhat*u(:);
        end
        
        % Sum performance measures
        perf(i) = perf(i) + sum(err.^2) + sum(ctr.^2);
        
        % Simulate one more step
        idx = Agents(1).pattern_idx;
        t   = sim.step();
    end

    % Deallocate network object. Required for SinrNetwork
    delete(Network)

    % Check whether the simulation converged prematurely
    conv(i) = t < simconf.Tf;

    % Update progress meter
    meter.notify(i);
end

%% Evaluate
% Check whether the simulations for each netconf converged
converged = all(conv, 3:4);

% Calculate mean performance for each netconf
h2 = sqrt(sum(perf, 3:4) / simconf.samples);
end

%% Helper functions
function err = formation_error(Agents)
    pos = [Agents.position];
    ref = [Agents.ref];
    err = norm(pos - mean(pos - ref,2) - ref, 'fro');
end
