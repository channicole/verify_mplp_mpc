% Author:   Nicole Chan
% Created:  11/10/17
% Description: This is the template user-defined input file
%
function funHand=inputParams()
    % Create function handles to the hybrid model components, MPC
    % formulation, and verification parameters
    funHand.flowEq = @flowModel;        % ODEs in each discrete mode (possibly nonlinear)
    funHand.modelJac = @modelJacobian;  % Jacobian in each mode
    funHand.unsafeStates = @unsafeSet;  % Sets of unsafe states
    funHand.initStates = @initialSet;   % Initial set of states for verification
    funHand.verifyParams = @verificationParams; % Time horizon, error terms, sim step-size
    funHand.MPCprob = @MPCprob;         % Linear, discrete-time plant model used to design MPC sol
end

function dX=flowModel(t,X,Fi,Gi,mode)
    % Example 4
    A = [1,1;0,1];
    B = [0;1];
    C = [0.5,0.5];
    dX=zeros(2,1);

    %% Closed-loop system
    ui = simTrace.Fi*X+simTrace.Gi;
    dX = eye(2)*X + Bsim*ui; % assumes we only implement the first optimal control step always

end

function Jac=modelJacobian(X,mode)
    Jac = [];
end

function filename = MPCprob()
    % The linear discrete-time approx. of the system must be given here
    % Give cost function, optimization horizon, and basic inf-norm
    % constraints on state/input/output vectors
    
    dim_x = 2;  % state vector length
    dim_u = 1;  % input vector length
    dim_y = 1;  % output vector length
    % TODO: read about choosing Ny and Nu << Ny: https://www.mathworks.com/help/mpc/ug/choosing-sample-time-and-horizons.html
%     Nc = 2;     % constraint horizon (number of future steps over which the constraints must hold)
    Nu = 2;     % control horizon (number of future nonzero control moves to solve for)
    Ny = 2;     % prediction horizon (number of future steps over which the cost is evaluated)

    % Cost function:
    % TODO: check dimensions of cost wrt dim_{} vars
    P = [];
    Q = [1,1;0,1];
    R = 0.8;
    
    % MPC constraints, format options:
    % % 1. element-wise min in first col, max in second col
    % % 2. [xmin,xmax] if using inf-norm for vector x
    % % 3. [] if no constraint
    % TODO: check dimensions wrt dim_{} vars
    xbnd = [-10*ones(dim_x,1),10*ones(dim_x,1)];
    ubnd = [-1,1];
    ybnd = [];
    
    % Linear model in state-space rep: x_dot = Ax + Bu, y = Cx
    A = [1,1;0,1];
    B = [0;1];
    C = [0.5,0.5];
    
    % Save to file and output the filename string
    filename = 'MPCmodel.mat';
    save(filename,'dim_x','dim_u','dim_y','Nu','Ny','P','Q','R','xbnd',...
        'ubnd','ybnd','A','B','C');
end

function Uset=unsafeSet()
    % Currently unsafe sets are represented by Polyhedron objects
    numProp = 4;                    % number of properties to check
    Uset(numProp,1) = Polyhedron(); % preallocate memory
    for i=1:numProp
        Uset(i) = Polyhedron();     % TODO: replace with actual unsafe sets
    end
end

function X0=initialSet()
    X0 = Polyhedron();
end

function [Thorizon,simStep,delta,maxPart,LipConst,epsilonConst]=verificationParams() 
    Thorizon = [0,100];     % verification time horizon
    simStep = 0.01;         % (fixed) simulation step size
    delta = [];             % sampling step size (same as what's used for discrete-time model and MPC period)
    maxPart = 10;           % maximum number of partitions taken per mode before returning an UNKNOWN result
    LipConst = 1;           % Lipschitz constant
    epsilonConst = 0;       % upper bound on simulation error
    
end