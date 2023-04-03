%---------------------------------------------------------------------------------------------------
% For Paper
% "Robust Performance Analysis for Time-Varying Multi-Agent Systems with Stochastic Packet Loss"
% by C. Hespe and H. Werner
% Copyright (c) Institute of Control Systems, Hamburg University of Technology. All rights reserved.
% Licensed under the GPLv3. See LICENSE in the project root for license information.
% Author(s): Christian Hespe
%---------------------------------------------------------------------------------------------------

clear

addpath('analysis', 'simulation', 'util')
addpath(genpath('mas-simulation/lib'))

%% Set up the problem
kappa = 0.1;
sysD  = ss(1, 1, 0, 0, 1);
sysC  = ss(-kappa, 0, 0, 0, 1);
sysP  = ss(0, 0, 1, 0, 1);

% Generate an appropriate set of graphs
graphs = {};
graphs{1} = circular_graph(20, 1, false);
graphs{2} = circular_graph(20, 2, false);
graphs{3} = circular_graph(20, 3, false);
graphs{4} = circular_graph(20, 4, false);
graphs{5} = circular_graph(20, 5, false);
graphs{6} = circular_graph(20, 6, false);
graphs{7} = circular_graph(20, 7, false);
N = height(graphs{1}.Nodes);

% Settings for the Monte-Carlo simulation
simconf = struct;
simconf.Tf        = 1e6;  % Maximum simulation time [s]
simconf.tol       = 1e-3; % Tolerance before stopping the simulation
simconf.samples   = 10;   % Number of samples per probability
simconf.R         = 0;    % Penalty on control signal
simconf.dT        = 1;    % Sampling time during simulation
simconf.positions = 1:N;  % Initial positions of the agents

% Additional properties of the communication links
range = Inf;  % Maximum range of communication
sym   = false; % If set to true, failure is always symmetric

%% Prepare synthesis and analysis steps
% Grid the probability axis. Higher density where larger changes are
% expected
p_swp = [0.0025:0.0025:0.0075, 0.01:0.005:0.03, 0.04:0.01:0.1, 0.12:0.02:0.3, 0.34:0.04:0.8, 0.82:0.02:1]';

% Prepare a set of switching sequences for the interconenction structure
netconf = struct('p', num2cell(p_swp), 'range', range, 'sym', sym);
pattern = [
    % Constant graph of each type
    1*ones(1,21);
    2*ones(1,21);
    3*ones(1,21);
    4*ones(1,21);
    5*ones(1,21);
    6*ones(1,21);
    7*ones(1,21);

    % Sequential switching
    repmat(1:7,1,3);

    % A few random switching sequences
    randi(length(graphs), 7, 21)
];

% Calculate eigenvalue bounds
lambda = [Inf, -Inf];
for graph = graphs
    l = eig(laplace_matrix(graph{:}));
    lambda(1) = min(lambda(1), l(2));
    lambda(2) = max(lambda(2), l(N));
end

%% Analyse the system
% Theoretical upper bound
h2_rb = zeros(size(p_swp));
tic
parfor i = 1:length(p_swp)
    h2_rb(i) = h2norm_decomposed_robust(sysD, sysC, sysP, N, lambda, p_swp(i));
end
disp(['Robust analysis completed in ' format_duration(toc)])

% Monte-Carlo analysis
tic
[h2_mc, converged] = mc_simulate(graphs, pattern, simconf, netconf);
disp(['Monte-Carlo analysis completed in ' format_duration(toc)])
converged = all(converged, 2);

% Filter out results that did not converge
h2_mc(~converged,:) = NaN;

% Determine min and max values
h2_min = min(h2_mc, [], 2);
h2_max = max(h2_mc, [], 2);

figure()
plot(p_swp, h2_rb, p_swp, h2_min, p_swp, h2_max)
hold on
patch([p_swp' fliplr(p_swp')], [h2_max' fliplr(h2_min')], 'g')
hold off
xlim([0,1])
ylim padded
xlabel('Transmission probability p')
ylabel('H_2 Performance')
legend('Robust', 'Monte-Carlo')

%% Export results
name = sprintf('performance_span_%d.csv', uint32(posixtime(datetime())));
tbl = table(p_swp, h2_rb, h2_min, h2_max);
writetable(tbl, name)
