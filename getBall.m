% Author:   Nicole Chan
% Created:  11/10/17
% Description: Takes a radius and outputs a ball (tbd which norm and which
% data struct) of that radius centered at 0.
%
function ball=getBall(dim,radius,center)
    if length(radius)==1
        radius = radius*ones(dim,1);
    end
    if nargin == 3 && length(center) == dim
        if isrow(center)
            center = center';
        end
        if isrow(radius)
            radius = radius';
        end
        ball = Polyhedron('lb', center-radius, 'ub', center+radius);
    elseif nargin == 2
        if isrow(radius)
            radius = radius';
        end
        ball = Polyhedron('lb', -radius, 'ub', radius);
    else
        error('Check input arguments.');
    end
%     % Compute all vertices of the box centered at origin (basically
%     % combination with replacement)
%     vals = [-radius,radius];
%     X = cell(1, dim);
%     [X{:}] = ndgrid(vals);
%     X = X(end : -1 : 1); 
%     vert = cat(dim+1, X{:});
%     vert = reshape(vert, [2^dim, dim]);
%     
%     % Create polyhedron object with the ball
%     ball = Polyhedron(vert);
end