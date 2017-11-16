% Author:   Nicole Chan
% Created:  11/10/17
% Description: Takes a Polyhedron object, overapproximates the polyhedron
% with a ball (either L1 or Linfty if using polyhedral data structs or L2
% if using ellipsoids), and returns the resulting ball's center and
% tightest radius achievable
%
% TODO: return center state as a row vector (column corresponds to number
% of intersecting regions)
function [center,radius,outPoly]=poly2ball(inPolyhedron,invariantRegions)
    if nargin==1
        % Currently uses the bounding box method to get inf-norm ball cover
        outPoly = inPolyhedron.outerApprox();

        % Compute center state from the ball and return largest radius
        center = zeros(1,inPolyhedron.Dim());
        radius = 0;
        for i=1:inPolyhedron.Dim()
            rad = max(outPoly.V(:,i))-min(outPoly.V(:,i))/2;
            center(i) = min(outPoly.V(:,i)) + rad;
            if rad > radius
                radius = rad;
            end
        end
    elseif nargin==2
        % Intersect inPolyhedron with invariantRegions 
        
    end
    
end