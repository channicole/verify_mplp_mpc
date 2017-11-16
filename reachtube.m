% Author:   Nicole Chan
% Created:  10/6/17
% Description: The 'reachtubeObj' object is a subclass of the doubly-linked 
% list node. It represents a polyhedral subset of the state-space that
% contains all the reachable states over the time horizon 'T'. More
% specifically, it contains all the reachsets R_i, for all i in T, where 
% R_i is the reachable states over t in [t_{i-1},t_{i}).
%
% Each reachtubeObj that gets generated has its own ID. An initial set
% of states may be partitioned to create sub-reachtubes, which are referred
% to as children of the original parent reachtubeObj object.
%
classdef reachtubeObj < dlnode
    properties
        Theta   % initial set of states
        x0      % initial state in Theta to simulate from
        xi      % nominal trajectory simulated
        MPCi    % affine control solution applied
        Reach   % reachsets (Reach(T(1))=Theta)
        T       % discrete time steps
        dia     % dim-dimensional vector of cover diameters
        ID      % array position if array is used to contain pointers
        parID   % parent node's ID
        childID % array of children IDs
    end
    
    methods
        function c = reachtubeObj(Theta,x0,xi,MPCi,Reach,T,dia,ID,parID) % constructor
            c = c@dlnode();
            if nargin == 0
                c.Theta = [];
                c.x0 = [];
                c.xi = [];
                c.MPCi = [];
                c.Reach = [];
                c.T = [];
                c.dia = [];
                c.ID = [];
                c.parID = [];
                c.childID = [];
            elseif nargin == 9
                c.Theta = Theta;
                c.x0 = x0;
                c.xi = xi;
                c.MPCi = MPCi;
                c.Reach = Reach;
                c.T = T;
                c.dia = dia;
                c.ID = ID;
                c.parID = parID;
                c.childID = [];
            else
                error('Incorrect number of arguments')
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
            dim = length(cov.dia);
            newDia = (cov.Yup(i,1:dim)-cov.Ylow(i,1:dim))./2;
            
            cov.x0 = cov.Y(i,:);
            cov.t0 = cov.T(i);
            cov.Y = cov.Y(i,:);
            cov.Yup = cov.Yup(i,:);
            cov.Ylow = cov.Ylow(i,:);
            cov.T = cov.T(i);
            cov.dia = newDia;
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
            dim = length(cov.dia);
            cov.Yup = max(cov.Yup,newCov.Yup(i,:));
            cov.Ylow = min(cov.Ylow,newCov.Ylow(i,:));
            cov.dia = (cov.Yup(1,1:dim)-cov.Ylow(1,1:dim))./2; % radius
            cov.T = max(cov.T,newCov.T(i));
            cov.t0 = cov.T(1);
            cov.x0 = cov.Ylow(1,1:dim)+cov.dia; % need to update support vars
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
%             cov.dia = vertcat(cov.dia,newCov.dia);
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
            cov.dia = newCov.dia;
            cov.children = [];
        end
    end % end methods
end % end classdef