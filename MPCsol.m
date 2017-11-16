% Name:     MPCsol_ARCaDiA.m
% Author:   Nicole Chan
% Date created: 10/3/17
% Description: Loads car models for different parameters, solves MPC
% solution, 

function [Pn,Fi,Gi,activeConstraints,Phard,details]=MPCsol(T_horizon,T_sampSize)
% Setup
if nargin == 1
    T_sampSize = 0.1;   % let sampling size be 0.1 sec
elseif nargin ~=2
    T_sampSize = 0.1;   % let sampling size be 0.1 sec
    T_horizon = 4;      % 4 sec
end

model = getCarModelForMPC(T_horizon,T_sampSize);

% Solve for MPC solution
[~,mplp_sol]=MPC_to_mpLP(model);

Pn = mplp_sol.xopt.Set;
% ACTIVECONSTRAINTS ONLY COMPUTED USING 'MPLP' SOLVER NOT 'PLCP'
try 
    activeConstraints = mplp_sol.mplpsol.activeConstraints;
catch
    activeConstraints = {};
end
Phard = mplp_sol.xopt.Domain;

% Optimizer
Fi = cell(1,mplp_sol.xopt.Num); 
Gi = cell(1,mplp_sol.xopt.Num);
if mplp_sol.xopt.Num>0
    optimizer = mplp_sol.xopt.Set.getFunction('primal');
    for i = 1:length(optimizer)
        Fi{i} = optimizer(i).F;
        Gi{i} = optimizer(i).g;
    end
end

% Cost
details.Ai = cell(1,mplp_sol.xopt.Num);
details.Bi = cell(1,mplp_sol.xopt.Num);
details.Ci = cell(1,mplp_sol.xopt.Num);
if mplp_sol.xopt.Num>0
    cost = mplp_sol.xopt.Set.getFunction('obj');
    for i = 1:length(cost)
        if isa(cost(i), 'QuadFunction')
            details.Ai{i} = cost(i).H;
        end
        details.Bi{i} = cost(i).F;
        details.Ci{i} = cost(i).g;
    end
end

end