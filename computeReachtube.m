% Author:   Nicole Chan
% Created:  10/9/17
% Description: Computes the reachtube of the given initial set of states
% over the given time horizon for a continuous (single-mode) system:
% 1. Computes  a 'center' state to simulate from
% 2. Simulate for length of time T
% 3. Compute discrepancy function and bloat the simulation trace
%
% TODO: make sure modelMPC is a struct with fields: Pn, F, G
function [Reachtube,safeflag]=computeReachtube(Theta,T,modelFun,modelMPC,modelJac,unsafeSet,epsilonConst,LipConst)
    %% Setup
    if nargin == 6
        epsilonConst = 0;
        LipConst = 0;
    elseif nargin == 7
        LipConst = 0;
    end

    %% Get ball approximation of the initial set
    [x0,deltaConst,ball0] = poly2ball(Theta,modelMPC.Pn);
    
    numRegions = length(deltaConst);            % Get number of CR's that intersect with Theta
    Reachtube = reachtubeObj();

    %% Compute a Reachtube for each discrete mode in Theta
    for j=1:numRegions
        %% Simulate
        % TODO: FIX ME.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Get control params F, G (u = F*x + G)
        CR = find(modelMPC.Pn.contains(x0(j,:)'));
        if ~isempty(CR)
            Fi = modelMPC.F{CR(1)}(Ny+Nu+1:Ny+Nu+dim_u,:);  % index CR in case there's multiple
            Gi = modelMPC.G{CR(1)}(Ny+Nu+1:Ny+Nu+dim_u);    % TODO: pass in these params!!
        else
            disp('Infeasible region reached');
            safeflag = 0;
            return
        end
        [sim.T,sim.X] = ode45(@(t,x) modelFun(t,x,Fi,Gi),T,x0); % TODO: check user-defined sim step size is set
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        [sim.T,sim.X] = simulateMPC();
        
        %% Compute discrepancy (ComputeLDF algorithm from ATVA15)
        b = zeros(length(sim.T),1);
        b(1) = deltaConst(j);
        beta = b;
        for i=2:length(sim.T)
            dt = sim.T(i)-sim.T(i-1);                   % length of timestep
            dia = (Delta + epsilonConst)*exp(LipConst*dt);% amount of coarse bloating
            S = Polyhedron([sim.X(i-1,:);sim.X(i,:)]);  % coarse overapprox of Reach over dt
            if dia ~= 0
                S = plus(S,getBall(dia));
            end
            [centerS,~,~] = poly2ball(S);               % get center state of S
            J = modelJac(centerS);                      % get Jacobian
            lambda = max(eig(J+J')/2);                  % 
            err = 0;                                    % ignore line 8 for now
            b(i) = lambda + err/2;                      % store local discrepancy
            Delta = (Delta + epsilonConst)*exp(b(i)*dt);% update error
            beta(i) = beta(i-1)*exp(b(i)*dt);           % missing from Alg 2
        end

        %% Bloat simulation trace
        u.Fi = Fi;
        u.Gi = Gi;
        rtube = reachtubeObj(ball0,x0,sim.X,u,[],sim.T,deltaConst(j),[],[]);
        rtube.Reach(1) = plus(sim.X(1,:),getBall(beta(1)));
        for i=2:length(sim.T)
            rtube.Reach(i) = plus(sim.X(i,:),getBall(beta(i)));
        end
        
        %% Check safety
        safeflag = checkSafety(rtube,unsafeSet);
        if ~safeflag
            return
        end
        
        %% Update Reachtube
        Reachtube.union(rtube);
    end
end