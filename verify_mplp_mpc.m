% Author:   Nicole Chan
% Created:  10/6/17
% Description: Implements C2E2's algorithm for a specific class of hybrid
% models, which use the (linear) MPC formulated using a mp-LP solver. 
% The model transitions periodically (as determined by sampling size in
% MPC) and the number of discrete modes correspond to the number of regions
% in the MPC solution.
% The controller is linear state feedback, but the plant dynamics may be
% nonlinear.
% For now, we only "refine covers" of the initial set, not subsequent
% reachsets. Furthermore, we do not search for counterexamples and decouple
% this function from safety certification.
function safe=verify_mplp_mpc(inputFile)
%% Setup
safe = -1;                          % -1: UNKNOWN, 0: UNSAFE, 1: SAFE

% Load user input from the file in the input arg
inParams = inputFile();             % instantiates function handles to each input entity needed from user
inFunNames = fieldnames(inParams);  % get function handle names from user input, stored in a cell array
if length(inFunNames) ~= 6          % TODO: update with inputParams
    error('Incorrect number of function handles. Follow template for the inputParams() method.');
else
    template = {'flowEq';'modelJac';'unsafeStates';'initStates';'verifyParams';'MPCprob'}; % TODO: update with inputParms
    if isempty(setdiff(inFunNames,template))
        error('Incorrect function handle names. Follow template for the inputParams() method.');
    end
end

% TODO: Check the parameters in each function for correctness

%% Compute the MPC solution
[~,OptSol] = MPC_to_mpLP(inParams.MPCprob);
[Pn,F,G,~,~,~] = MPCsol(OptSol);

%% Safety Verification
tic;
% Generate reachtubes using DFS
covers = stack(Theta);
while (~covers.isempty)
    % Get cover at the top of stack (don't pop until it's been checked)
    currCov = covers.getLast();
    
    % If multiple modes are valid in Theta, partition into separate covers
    % and add to stack (don't forget to update child/parIDs
    if currCov.Theta
    
    % Initialize Theta to top of stack
    
    % If length of T > Tsamp, then duplicate cover, change T(1) and push to
    % stack and update current Theta
    
    % Call computeReachtube for T(1):T(1)+Tsamp (if T(2)-T(1) > Tsamp)
    
    % Push last reachset returned and adjust T(1):=T(1)+Tsamp (this becomes
    % child of the current set) (if T(2)-T(1) > Tsamp)
    
    % Else: find parent, copy reachtube into parent (take union with 
    % pre-exiting data if applicable), and pop
    
    
end % end verification no covers left to check
toc;

end