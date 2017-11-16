% Name:     simulateMPC.m
% Author:   Nicole Chan
% Date created: 11/15/17 (modified from ARCaDiA_MPC repo)
% Description: Loads DT-model, computes MPC solution offline, and simulates
% the CT-system implementing the MPC. Preparation for verification.

% 'controlParamFile' is what will be passed to the MPC solver
% 'simParamFile' contains initial condition, model (possibly different
% from the linear one used for MPC), time horizon, and time step size for 
% simulation and will also contain constraints for verification

% NOTE: currenly only simulating/recording x, not y even if the output IS
% specified. Don't forget to update when there's a need.
% NOTE: currently assumes sampling time delta=1, so delta never shows up.
% Fix this by adding delta as a parameter to the discrete time model in
% 'controlParamFile'

function [sim_T,sim_X]=simulateMPC(controlParamFile,simParamFile)
%% Setup
    safeFlag = 1;
    options = odeset('RelTol',1e-8,'AbsTol',1e-8); % options for ode45
    err = warning; % save current warning settings
    warning('error','MATLAB:load:variableNotFound'); % set this warning to throw error

    if nargin ~= 2
        error('Incorrect input arguments.');
    else
        try
            load(controlParamFile,'A','B','C','dim_u','dim_x','dim_y',...
                'Nc','Nu','Ny','P','Q','R','ubnd','xbnd','ybnd');
        catch
            error('Cannot load MPC problem parameters, parameters may be incomplete.');
        end
        try
            load(simParamFile,'Asim','Bsim','Csim','Tf','Tstep','x0');
        catch
            error('Cannot load simulation parameters, parameters may be incomplete.');
        end
    end

    % simParamFile must include a function handle for the ODEs modeling the
    % system if nonlinear. Otherwise, linear models can use a function handle
    % directly or be in state-space model form
    try
        load(simParamFile,'modelFunHandle');
        f = modelFunHandle;
    catch
        f = @modelFun;
    end

    % Create ODE function if a state-space model is specified
    function dX=modelFun(t,X)
        dX=zeros(dim_x,1);

        %% Closed-loop system
        ui = simTrace.F{i}*X+simTrace.G{i};
        dX = Asim*X + Bsim*ui; % assumes we only implement the first optimal control step always
    end

%% Compute MPC controller
    [Pn,Fi,Gi,~,~,~]=MPCsol_mpLP_solve(controlParamFile,'optSol1.mat');
    if isempty(Pn)
        error('No feasible solution was found.');
    end

%% Simulate
    simTrace.t = (0:Tstep:Tf)';
    simTrace.x = zeros(Tf/Tstep+1,dim_x);
    simTrace.x(1,:) = x0';
    simTrace.u = zeros(Tf/Tstep+1,1);
    % Control law for each update is u(t) = F*x(t)+G
    simTrace.F = cell(Tf,1);
    simTrace.G = cell(Tf,1);

    for i=1:Tf
        j = (i-1)/Tstep+1;
        % Update control law
        x0 = simTrace.x(j,:)';
        CR = find(Pn.contains(x0));
        if ~isempty(CR)
            simTrace.F{i} = Fi{CR(1)}(Ny+Nu+1:Ny+Nu+dim_u,:); % index CR in case there's multiple
            simTrace.G{i} = Gi{CR(1)}(Ny+Nu+1:Ny+Nu+dim_u);
        else
            error('Infeasible region reached');
        end

        % Simulate by solving ODE
        [~,Xout] = ode45(f,simTrace.t(j:j+1/Tstep),x0,options);

        % Update trace variables
        simTrace.x(j:j+1/Tstep,:) = Xout;
%         Uout = Xout*simTrace.F{i}'+repmat(simTrace.G{i}',size(Xout,1),1);
%         simTrace.u(j:j+1/Tstep,:) = Uout(:,Ny+Nu+1:Ny+Nu+dim_u);

    end

%% Closing statements
warning(err); % return warning settings to default

end