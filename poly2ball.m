% Author:   Nicole Chan
% Created:  11/10/17
% Description: Takes a Polyhedron object, overapproximates the polyhedron
% with a ball (either L1 or Linfty if using polyhedral data structs or L2
% if using ellipsoids), and returns the resulting ball's center and
% tightest radius achievable
%
% 'invariantRegions' (or 'Pn' from MPC solution) is a Polyhedron-array
% 'inPolyhedron' must be a single Polyhedron (not array)
%
function [center,radius,outBall,outPoly]=poly2ball(inPolyhedron,invariantRegions)
    % Initialize output variables
    center = [];
    radius = [];
    outBall = [];
    outPoly = inPolyhedron;
    
    % If we just want the ball to cover the input polyhedron:
    if nargin==1
        % Currently uses the bounding box method to get inf-norm ball cover
        outBall = inPolyhedron.outerApprox();

        % Compute center state from the ball and return largest radius
        center = zeros(inPolyhedron.Dim(),1);
        radius = zeros(inPolyhedron.Dim(),1);
        for i=1:inPolyhedron.Dim()
            rad = (max(outBall.V(:,i))-min(outBall.V(:,i)))/2;
            center(i) = min(outBall.V(:,i)) + rad;
            radius(i) = rad;
        end
        center = center';
    % If we want to first check for intersection with different modes and
    % return corresponding ball-covers
    elseif nargin==2 && length(inPolyhedron)==1
        % Intersect inPolyhedron with invariantRegions 
        intRegions = invariantRegions & inPolyhedron;
        
        % Return array of polyhedra (non-empty intersections)
        outPoly = intRegions(intRegions.isFullDim());
        
        % Return array of ball-approximations
        outBall = outPoly.outerApprox();
        
        % Return array of center states and corresponding largest radii
        center = zeros(length(outPoly),inPolyhedron.Dim());
        radius = zeros(length(outPoly),inPolyhedron.Dim());
        for j=1:length(outPoly)
            for i=1:inPolyhedron.Dim()
                rad = (max(outBall(j).V(:,i))-min(outBall(j).V(:,i)))/2;
                center(j,i) = min(outBall(j).V(:,i)) + rad;
                radius(j,i) = rad;
            end
        end
    end  
end