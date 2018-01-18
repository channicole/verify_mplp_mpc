% Author:   Nicole Chan
% Created:  10/9/17
% Description: Computes the reachtube of the given initial set of states
% over the given time horizon for a continuous (single-mode) system:
% 1. Computes  a 'center' state to simulate from
% 2. Simulate for length of time T (T should be a vector of all simulation
% time steps)
% 3. Compute discrepancy function and bloat the simulation trace
%
% TODO: make sure modelMPC is a struct with fields: Pn, F, G
function [Reachtube,safeflag]=computePost(Theta,T,modelFun,modelMPC,modelJac,unsafeSet,epsilonConst,LipConst)
    %% Setup
    if nargin == 6
        epsilonConst = 0;
        LipConst = 0;
    elseif nargin == 7
        LipConst = 0;
    end

    % Get ball approximation of the initial set
    [x0,deltaConst,ball0] = poly2ball(Theta,modelMPC.Pn);
    % % % Sanity check % % %
    if size(deltaConst,1)~=1; error('In computePost, the initial set should only correspond with one invariant region.'); end;
    % % % % % % % % % % % % 
    Reachtube = reachtubeObj(Theta,x0,x0',[],ball0,T,deltaConst);

    %% Simulate
    % Get control params Fi, Gi (u = Fi*x + Gi)
    CR = find(modelMPC.Pn.contains(x0'));
    if ~isempty(CR)
        Reachtube.MPCi.Fi = modelMPC.Fi{CR(1)};                     % index CR in case there's multiple
        Reachtube.MPCi.Gi = modelMPC.Gi{CR(1)};
    else
        disp('Infeasible region reached');
        safeflag = 0;
        return
    end
    [sim.T,sim.X] = ode45(@(t,x) modelFun(t,x,Reachtube.MPCi.Fi,Reachtube.MPCi.Gi),T,x0); 

    %% Compute discrepancy (ComputeLDF algorithm from ATVA15)
    x_dim = length(x0);
    b = zeros(length(sim.T),1);
    b(1) = max(deltaConst);
    beta = b;
    Delta = b(1);
    for i=2:length(sim.T)
        dt = sim.T(i)-sim.T(i-1);                   % length of timestep
        dia = (Delta + epsilonConst)*exp(LipConst*dt);% amount of coarse bloating
        S = Polyhedron([sim.X(i-1,:);sim.X(i,:)]);  % coarse overapprox of Reach over dt
        if dia ~= 0
            S = plus(S,getBall(x_dim,dia));
        end
        [centerS,~,~] = poly2ball(S);               % get center state of S
        J = modelJac(centerS,Reachtube.MPCi.Fi,Reachtube.MPCi.Gi); % get Jacobian
        lambda = max(eig(J+J')/2);                  % 
        err = 0;                                    % ignore line 8 for now
        b(i) = lambda + err/2;                      % store local discrepancy
        Delta = (Delta + epsilonConst)*exp(b(i)*dt);% update error
        beta(i) = beta(i-1)*exp(b(i)*dt);           % missing from Alg 2
    end

    %% Bloat simulation trace
    rtube(length(sim.T),1) = Polyhedron();
    for i=1:length(sim.T)
        rtube(i) = getBall(x_dim,beta(i),sim.X(i,:));
    end

    %% Check safety
    safeflag = checkSafety(rtube,unsafeSet);
    if ~safeflag
        return
    end

    %% Update Reachtube parameters with computed reachtube data
    Reachtube.updateReach(rtube);
    Reachtube.MPCi = repmat(Reachtube.MPCi,length(sim.T),1);
    Reachtube.xi = sim.X;
end