clear

addpath('analysis', 'util')

%% Set up the problem
kappa = 0.1;

% Consensus Example
conD = ss(1, 1, 0, 0, 1);
conC = ss(-kappa, 0, 0, 0, 1);
conP = ss(0, 0, 1, 0, 1);

% Input/Output-swapped system
adjD = ss(1, 0, 1, 0, 1);
adjC = ss(-kappa, 0, 0, 0, 1);
adjP = ss(0, 1, 0, 0, 1);

N   = 20;   % Number of agents
p   = 0.5;  % Transmission probability
slc = 5;    % lambda_min for which the slice of the function is taken

points     = 200;   % Number of gridpoints on each dimension
lambda_min = 0.08;  % Lower bound on the lambda grid
lambda_max = N;     % Upper bound on the lambda grid

% Generate grid to evaluate the function on
x = unique(sort([linspace(lambda_min, lambda_max, points), slc]));
y = linspace(lambda_min, lambda_max, points);

[X, Y] = meshgrid(x, y);
Z_con  = NaN(size(X));
Z_adj  = NaN(size(X));

%% Calculate upper bound on H2 norm for each grid point
parpool([8, 32]);

% Create progress meter
meter = ParforProgressMeter(numel(X), 0.02);
meter.start()

tic
parfor i = 1:numel(X)
    % Skip if lambda_min > lambda_max
    if X(i) <= Y(i)
        Z_con(i) = h2norm_decomposed_robust(conD, conC, conP, N, [X(i), Y(i)], p);
        Z_adj(i) = h2norm_decomposed_robust(adjD, adjC, adjP, N, [X(i), Y(i)], p);
    end
        
    meter.notify(i);
end
disp(['Parameter sweep completed in ' format_duration(toc)])

% Shutdown parpool
delete(gcp('nocreate'));

%% Generate contour lines
% The contour lines are generated on a log scale such that they are not all
% squished together in the top left corner. The labelling will be reverted
% to linear later.

figure()
data_con = contour(X, Y, Z_con, 20:10:90)';
title('Consensus Problem')
xlabel('\lambda_{min}')
ylabel('\lambda_{max}')

figure()
data_adj = contour(X, Y, log10(Z_adj), 1.2:0.2:2.6)';
title('Adjoint System')
xlabel('\lambda_{min}')
ylabel('\lambda_{max}')

%% Generate vertical slice of the data
idx  = find(x == slc);
mask = ~isnan(Z_con(:,idx));
y_slc = Y(mask, idx);
h2_slc = Z_con(mask, idx);

figure()
plot(y_slc, h2_slc)
xlabel('\lambda_{max}')
ylabel('Bound on the H_2 norm')
title(['Slice of the Performance field for \lambda_{min} = ' num2str(slc)])

%% Export the data
name = sprintf('eigenvalue_contour_%d.csv', uint32(posixtime(datetime())));
writematrix(data_con, name);

name = sprintf('eigenvalue_contour_adj_%d.csv', uint32(posixtime(datetime())));
writematrix(data_adj, name);

name = sprintf('eigenvalue_slice_%d.csv', uint32(posixtime(datetime())));
tbl = table(y_slc, h2_slc);
writetable(tbl, name)
