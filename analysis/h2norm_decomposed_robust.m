function [H2, Q, solver_stats] = h2norm_decomposed_robust(sysD, sysC, sysP, N, lambda, p)
%H2NORM_DECOMPOSED_ROBUST Calculate an upper bound on the H2-norm of a
%time-varying decomposable jump system in a scalable manner
%   This function implements Theorem 7 of the paper the code belongs to.
%
%   Arguments:
%       sysD   -> Decoupled part of the system
%       sysC   -> Stochastically coupled part of the system
%       sysP   -> Deterministically coupled part of the system
%       N      -> Number of agents in the multi-agent system
%       lambda -> Bounds on the Laplacian spectrum
%       p      -> Probability of a successful package transmission
%   Returns:
%       H2           -> Upper bound on the H2-norm of the system
%       Q            -> Storage function matrix of the solution
%       solver_stats -> Timing information about the algorithm

% Setup the SDP solver
% The offset is used to convert strict LMI into nonstrict ones by pushing
% the solution away from the boundary.
offset = 1e-8;
opts   = sdpsettings('verbose', 0);

% Allocate storage for time measurements
solver_stats        = struct;
solver_stats.prep   = NaN;
solver_stats.solver = NaN;
solver_stats.total  = NaN;
solver_stats.yalmip = NaN;

% Initialize return values with reasonable defaults
H2    = inf;
timer = tic;

%% Prepare the system description
[Ad, Bd, Cd, Dd] = ssdata(sysD);
[Ac, Bc, Cc, Dc] = ssdata(sysC);
[Ap, Bp, Cp, Dp] = ssdata(sysP);

nx = size(Ac,1);
nw = size(Bc,2);

%% Solve the SDP from Theorem 2
Q = sdpvar(nx, nx);
Z = sdpvar(nw, nw);

% Initialize cost and constraints
Constraints = Q >= offset * eye(nx);
cost        = trace(Z);

% The eigenvalues of L0 are sorted by Matlab. By iteration only over 2:N,
% we ignore lambda_1 = 0 and thus the uncontrollable and marginally stable
% modal subsystem.
for i = 1:2
    li  = lambda(i);
    lit = p*(1-p)*li;
    
    Acl = Ad + p*li*Ac + li*Ap;
    Bcl = Bd + p*li*Bc + li*Bp;
    Ccl = Cd + p*li*Cc + li*Cp;
    Dcl = Dd + p*li*Dc + li*Dp;
    
    LMI = Acl'*Q*Acl + Ccl'*Ccl - Q + 2*lit * (Ac'*Q*Ac + Cc'*Cc);
    TRC = Bcl'*Q*Bcl + Dcl'*Dcl - Z + 2*lit * (Bc'*Q*Bc + Dc'*Dc);
    
    % For certain LMIs that are obviously infeasible, Yalmip will refuse to
    % construct the SDP and issue an error. This try catch will convert
    % that error into a warning so that we can successfully finish running
    % the remainder of the script.
    try
        Constraints = [ Constraints                     ,...
                        LMI <= -offset * eye(size(LMI)) ,...
                        TRC <= -offset * eye(size(TRC)) ];
    catch ME
        warning(ME.message)
        Q = [];
        solver_stats.total = toc(timer);
        return
    end
end

solver_stats.prep = toc(timer);
sol = optimize(Constraints, cost, opts);

if sol.problem ~= 0
    warning('YALMIP return an error: %s', sol.info)
    Q = [];
else
    % We calculate gamma^2 with the LMI constraints, so we need to take the
    % square root here.
    H2 = sqrt((N-1)*value(cost));
    
    Q  = value(Q);
    solver_stats.yalmip = sol.yalmiptime;
    solver_stats.solver = sol.solvertime;
end

solver_stats.total = toc(timer);
end
