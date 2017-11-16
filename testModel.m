% Author:   Nicole Chan
% Created:  11/6/17
% Description: An example ODE model file
%
function dX=testModel(t,X,Fi,Gi)
dX=zeros(2,1);

%% Closed-loop system
ui = simTrace.Fi*X+simTrace.Gi;
dX = eye(2)*X + Bsim*ui; % assumes we only implement the first optimal control step always

end