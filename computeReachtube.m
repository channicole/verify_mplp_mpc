% Author:   Nicole Chan
% Created:  10/9/17
% Description: Computes the reachtube of the given initial set of states
% over the given time horizon for a continuous (single-mode) system:
% 1. Computes  a 'center' state to simulate from
% 2. Simulate for length of time T
% 3. Compute discrepancy function and bloat the simulation trace
%
function [reachtube,safeflag]=computeReachtube(Theta,modelFun,T)

% Simulate
[sim.T,sim.X] = ode45(@(t,x) modelFun(t,x,Fi,Gi),T,x0);

end