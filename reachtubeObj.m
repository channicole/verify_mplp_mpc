% Author:   Nicole Chan
% Created:  10/6/17
% Description: The 'reachtubeObj' object represents a polyhedral subset of 
% the state-space that contains all the reachable states over the time 
% horizon 'T'.
% It contains all the reachsets R_i, for all i in T, where 
% R_i is the reachable states over t in [t_{i-1},t_{i}).
%
classdef reachtubeObj < handle
    properties (Access = public)
        Theta   % Polyhedron: initial set of states
%         x0      % n-length vector: initial state in Theta to simulate from
        xi      % k-length array of n-length vectors: nominal trajectory simulated
        MPCi    % k-length array of MPC-objects: affine control solution applied
%         Reach   % k-length array of Polyhedron: reachsets (Reach(T(1))=Theta)
        T       % k-length vector: discrete time steps
%         rad     % n-length vector: radius of the ball approximation of Theta
    end
    % The following members represent the hyperrectangular approximation of
    % Theta. The Polyhedron representation is stored in Reach(1,:). Its
    % center state and radius is stored in x0 and rad, respectively.
    properties (Access = protected)
        x0      % n-length vector: initial state in Theta to simulate from
        Reach   % k-length array of Polyhedron: reachsets (Reach(T(1))=Theta)
        rad     % n-length vector: radius of the ball approximation of Theta
    end
    
    methods
        function reach = reachtubeObj(Theta,x0,xi,MPCi,Reach,T,rad) % constructor
            if nargin == 7
                reach.Theta = Theta;
                reach.x0 = x0;
                reach.xi = xi;
                reach.MPCi = MPCi;
                reach.Reach = Reach;
                reach.T = T;
                reach.rad = rad;
            else
                error('Incorrect number of arguments')
            end
        end
        
        function reach = updateReach(reach,newReach)
            if size(newReach,1) == 1
                reach.Reach = newReach.outerApprox();
                reach.rad = transpose(max(reach.Reach)-min(reach.Reach));
            elseif ~isempty(newReach)
                % % NOTE: currently assumes each Polyhedron-element of
                % newReach is already approximated by a ball/box
                reach.Reach = newReach;
                reach.rad = transpose(max(reach.Reach(1,:))-min(reach.Reach(1,:)));
            else
                warning('Input is an empty variable, so member not updated.')
            end
        end
        
        function reach = copyObject(inReach,reach)
            C = metaclass(inReach);
            P = C.Properties;
            for k=1:length(P)
                if ~P{k}.Dependent
                    reach.(P{k}.Name) = input.(P{k}.Name);
                end
            end
        end
        
        function c = union(c,inReach)
            if nargin < 2
                error('Incorrect number of arguments')
            end
            [~,ind1,ind2] = intersect(c.T,inReach.T);
            for i=1:length(ind1)
                c.Reach(ind1(i)).Union(inReach.Reach(ind2(i)));
            end
            ind2 = setdiff((1:length(inReach.T))',ind2);
            c.Reach(end+1:end+length(ind)) = inReach.Reach(ind2);
        end
        function cov = reduceCover(cov,i)
            % Creates new cover from the reachset at cov(i) 
            if nargin < 2
                error('Incorrect number of arguments')
            end
            dim = length(cov.rad);
            newDia = (cov.Yup(i,1:dim)-cov.Ylow(i,1:dim))./2;
            
            cov.x0 = cov.Y(i,:);
            cov.t0 = cov.T(i);
            cov.Y = cov.Y(i,:);
            cov.Yup = cov.Yup(i,:);
            cov.Ylow = cov.Ylow(i,:);
            cov.T = cov.T(i);
            cov.rad = newDia;
        end
        function cov = coverUnion(cov,newCov,i)
            % Updates cov to be union with newCov(i)
            % Used with single reachset cov for nextRegions
            % Used with multiple reachsets and vector i for stitching
            % together a reachtubeObj in computeReachtube
            if nargin < 3
                error('Incorrect number of arguments')
            end
            if length(i) ~= length(cov.T) || length(i) > length(newCov.T)
                error('Number of reachsets in cover does not match number of indices specified')
            end
            loc = cov.x0(end); % added to track nondet trans to passive
            dim = length(cov.rad);
            cov.Yup = max(cov.Yup,newCov.Yup(i,:));
            cov.Ylow = min(cov.Ylow,newCov.Ylow(i,:));
            cov.rad = (cov.Yup(1,1:dim)-cov.Ylow(1,1:dim))./2; % radius
            cov.T = max(cov.T,newCov.T(i));
            cov.t0 = cov.T(1);
            cov.x0 = cov.Ylow(1,1:dim)+cov.rad; % need to update support vars
            cov.x0 = ARPOD_update(cov.x0,cov.t0,loc); % help differentiate passive
            cov.Y = cov.x0; % need to update support vars
        end
        function addMode(cov,newCov) % newCov should occur sequentially after cov
            if nargin < 2
                error('Incorrect number of arguments')
            end
            % % MAY WANT TO ADD CHECK TO SEE THAT TIME IS UNIQUE AND
            % SEQUENTIAL
            cov.T = vertcat(cov.T,newCov.T);
            cov.Y = vertcat(cov.Y,newCov.Y);
            cov.Yup = vertcat(cov.Yup,newCov.Yup);
            cov.Ylow = vertcat(cov.Ylow,newCov.Ylow);
%             cov.rad = vertcat(cov.rad,newCov.rad);
        end
        function addChildren(cov,children)
            cov.children = children;
        end
        function numnodes = listlength(startNode)
            numnodes = 0;
            if isempty(startNode)
                return
            else
                head = startNode;
                numnodes = 1;
                while ~isempty(head.Next)
                    numnodes = numnodes+1;
                    head = head.Next;
                end
            end
        end
        function leaf = isleaf(cov)
            leaf = zeros(1,length(cov));
            % Check no children
            for i=1:length(cov)
                if isempty(cov(i).children)
                    if length(cov(i).T) <= 1
                        return
                    end
                    leaf(i) = 1;
                end
            end
        end
        function cov = copyCover(cov,newCov)
            % Copies newCov properties into cov, but maintains linked list
            % properties and ID/parent
            cov.x0 = newCov.x0;
            cov.t0 = newCov.t0;
            cov.Y = newCov.Y;
            cov.Yup = newCov.Yup;
            cov.Ylow = newCov.Ylow;
            cov.T = newCov.T;
            cov.rad = newCov.rad;
            cov.children = [];
        end
    end % end methods
end % end classdef