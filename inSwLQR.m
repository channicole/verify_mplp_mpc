% Author:   Nicole Chan
% Created:  1/13/18
% Description: Input file for (non-MPC) original two-stage LQR satellite
% example.
%
% NOTE: the unsafe set is actually specified by the safe set
%
function funHand=inSwLQR()
    % Create function handles to the hybrid model components, MPC
    % formulation, and verification parameters
    funHand.flowEq = @flowModel;        % ODEs in each discrete mode (possibly nonlinear)
    funHand.modelJac = @modelJacobian;  % Jacobian in each mode
    funHand.unsafeStates = unsafeSet;   % Sets of unsafe states
    funHand.initStates = initialSet;   % Initial set of states for verification
    funHand.vPar = verificationParams;  % Time horizon, error terms, sim step-size
    funHand.MPCprob = @MPCprob;         % Linear, discrete-time plant model used to design MPC sol
    funHand.MPCsol = swLQR;
end

function dX=flowModel(t,X,Fi,Gi)
    % Satellite plant dynamics
    n = 7.3023e-05*60;  % [1/min]
    m = 500;            % [kg]
    A = [0,0,1,0; 0,0,0,1;3*n^2,0,0,2*n;0,0,-2*n,0];
    B = [0,0;0,0;1/m,0;0,1/m];
    
    % Closed-loop system
    ui = Fi*X;          % Fi corresponds with Klqr from old code; Gi should be zeros, just ignore
    dX = A*X+B*ui;
end

function Jac=modelJacobian(X,Fi,Gi)
    % This example is LTI so no dependence on input argument X
    % Also the A-matrix from closed-loop dynamics becomes the Jacobian
    n = 7.3023e-05*60;  % [1/min]
    m = 500;            % [kg]
    A = [0,0,1,0; 0,0,0,1;3*n^2,0,0,2*n;0,0,-2*n,0];
    B = [0,0;0,0;1/m,0;0,1/m];
    Jac = A+B*Fi;
end

function filename = MPCprob()
    filename = [];
end

function sol=swLQR() 
    X_BND = 10000;  % [m] some large value of positions for invariant region
    X_GUARD = 100;  % [m] the magnitude of position the transition guard
    A_hex = [1,0,0,0;  1,1,0,0;...
             0,1,0,0; -1,1,0,0;...
            -1,0,0,0;-1,-1,0,0;...
             0,-1,0,0;1,-1,0,0];
    b_hex = X_GUARD.*[1;sqrt(2);1;sqrt(2);1;sqrt(2);1;sqrt(2)];
    
    R2 = Polyhedron('A',A_hex,'b',b_hex);
    R1 = Polyhedron('A',[1,0,0,0;0,1,0,0;-1,0,0,0;0,-1,0,0],'b',X_BND.*ones(4,1))-R2;
    % [~,F1,F2]=getlqr();
    F1 = [28.8286776769430,-0.100479948259883,1449.97541985328,-0.00462447231887482;...
        0.0870156786852279,33.2561992450513,-0.00462447231887482,1451.50134643428]; 
    F2 = [288.028766271474,-0.131243039715836,9614.98979543236,3.41199965400404e-07;...
        0.131243040368934,287.999970095943,3.41199965400404e-07,9614.98829796995];
    G1 = zeros(2,1);
    G2 = zeros(2,1);
    sol.Pn = [R1;R2];
    sol.Fi = {-F1,-F2};
    sol.Gi = {G1,G2};
end

function Uset=unsafeSet()
    % Currently unsafe sets are represented by Polyhedron objects
    numProp = 3;                            % number of properties to check
    Uset.region(numProp,2) = Polyhedron();  % preallocate memory
    Uset.safe = ones(3,1);       % 0: the region defines Unsafe, 1: the region defines Safe
    
    %% Max thrust
    U_BND = 36000;  % [kg*m/min^2] magnitude of safe thrust
    X_BND = 10000;  % [m] some large value of positions for invariant region
    X_GUARD = 100;  % [m] the magnitude of position the transition guard
    A_hex = [1,0,0,0;  1,1,0,0;...
             0,1,0,0; -1,1,0,0;...
            -1,0,0,0;-1,-1,0,0;...
             0,-1,0,0;1,-1,0,0];
    b_hex = X_GUARD.*[1;sqrt(2);1;sqrt(2);1;sqrt(2);1;sqrt(2)];
    
    R2 = Polyhedron('A',A_hex,'b',b_hex);
    R1 = Polyhedron('A',[1,0,0,0;0,1,0,0;-1,0,0,0;0,-1,0,0],'b',X_BND.*ones(4,1))-R2;
    % [~,F1,F2]=getlqr();
    F1 = [28.8286776769430,-0.100479948259883,1449.97541985328,-0.00462447231887482;...
        0.0870156786852279,33.2561992450513,-0.00462447231887482,1451.50134643428]; 
    F2 = [288.028766271474,-0.131243039715836,9614.98979543236,3.41199965400404e-07;...
        0.131243040368934,287.999970095943,3.41199965400404e-07,9614.98829796995];
    Uset.region(1,1) = (R1 & Polyhedron('A',[F1;-F1],'b',U_BND.*ones(4,1)))... 
        + (R2 & Polyhedron('A',[F2;-F2],'b',U_BND.*ones(4,1)));
    
    %% Total velocity in Phase 3
    V_BND = 3;      % [m/min] max magnitude of velocity to be safe
    A_hex = [0,0,1,0;  0,0,1,1;...
             0,0,0,1; 0,0,-1,1;...
             0,0,-1,0;0,0,-1,-1;...
             0,0,0,-1;0,0,1,-1];
    b_hex = V_BND.*[1;sqrt(2);1;sqrt(2);1;sqrt(2);1;sqrt(2)];
    Uset.region(2,1) = R2 & Polyhedron('A',A_hex,'b',b_hex);
    Uset.region(2,2) = R1;
    
    %% LOS
    A_tri = [-1,0,0,0;1/sqrt(3),1,0,0;1/sqrt(3),-1,0,0];
    b_tri = [100;0;0];
    Uset.region(3,1) = R2 & Polyhedron('A',A_tri,'b',b_tri);
    Uset.region(3,2) = R1;
end

function X0=initialSet()
    x0 = [-900,-400,0,0];
    rad = [5,5,0.1,0.1];
    X0 = getBall(4,rad,x0);
end

function vPar=verificationParams() 
    vPar.Thorizon = [0,270];     % verification time horizon
    vPar.simStep = 0.01;         % (fixed) simulation step size
    vPar.deltaStep = 10;         % sampling step size (same as what's used for discrete-time model and MPC period)
    vPar.maxPart = 4;            % maximum number of partitions taken per mode before returning an UNKNOWN result
    vPar.LipConst = 1;           % Lipschitz constant
    vPar.epsilonConst = 0;       % upper bound on simulation error
    
end