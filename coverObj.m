% Author:   Nicole Chan
% Created:  1/9/18
% Description: The 'coverObj' object is a subclass of 'reachtubeObj' that
% includes an extra parameter to indicate whether the entire reachtube has
% been computed and if so, if it is safe. The initial set of states for a
% 'coverObj' should always correspond to a partition of the user-defined
% inital set of states for verification purposes. The reachtube for a
% 'coverObj' should also cover the entire verification time horizon
% specified by the user. 
% In this way, 'coverObj' is distinguished from 'reachtubeObj', which is a 
% data structure that may be used for storing reachsets/tubes that only 
% cover some interval of the time horizon and use initial states that 
% result from subcomputations rather than the user-defined initial set
% Theta.
% 
% safeFlag: (1) reachtube computed and is safe; (0) reachtube computed and
% is not safe, so partitioning is needed; (-1) reachtube is not computed
% yet; (-2) exit verification process because the max number of partitions
% is reached without safe result
%
% NOTE: as of this iteration, we do not search for unsafe traces. 
%
classdef coverObj < reachtubeObj
    properties
        safeFlag   % indicates whether a reachtube has been computed yet and if it is safe
    end
    
    methods
        % Constructor
        function cover = coverObj(Theta,x0,xi,MPCi,Reach,T,rad)
            if nargin ~= 7
                Theta = [];
                x0 = [];
                xi = [];
                MPCi = [];
                Reach = [];
                T = [];
                rad = [];
            end
            cover = cover@reachtubeObj(Theta,x0,xi,MPCi,Reach,T,rad);
            cover.safeFlag = safetyEnum.NeedReach;
        end
        
        % Set safeFlag
        function cover = setSafeFlag(cover,flag)
            if nargin ~= 2
                error('Incorrect number of arguments')
            end
            %%%% TODO: CHECK IF REACH COMPUTED, CHECK IF ACTUALLY SAFE, UPDATE %%%%
            if cover.safeFlag ~= 1
                cover.safeFlag = flag;
            end
        end
        
        % Update the reachtube parameters
        function cover = updateReach(cover,inReach)
            if nargin ~= 2
                error('Incorrect number of arguments')
            end
            
            cover.x0 = inReach.x0;
            cover.xi = inReach.xi;
            cover.MPCi = inReach.MPCi;
            cover.Reach = inReach.Reach;
            cover.T = inReach.T;
            cover.rad = inReach.rad;
        end
        
        % Partition the current cover and returns the resulting 2^n covers
        % in a cell array
        function outCovers = partitionCover(inCover) 
            if nargin < 1
                error('Incorrect number of arguments')
            end
            
            n = length(inCover.x0);
            outCovers = cell(2^n,1);
            
            if length(inCover.x0)==n && length(inCover.rad)==n 
                x0 = inCover.x0;
                rad = inCover.rad;
            else
                % Get ball approximation of the initial set
                [x0,rad] = poly2ball(inCover.Theta);
            end

            if length(rad)~=n
                error('Yo you got dimension issues with member variable rad');
            end
            
            vals = [-0.5,0.5];
            X = cell(1, n);
            [X{:}] = ndgrid(vals);
            X = X(end : -1 : 1); 
            xcen = cat(n+1, X{:});
            xcen = reshape(xcen, [2^n, n]);
            xcen = xcen*diag(rad) + repmat(x0',2^n,1);
          
            rad = rad/2;
            for i=1:2^n
                % Get hyperbox that approximates the i-th partition of inCover.Theta
                boxTheta = Polyhedron('lb',xcen(i,:)'-rad,'ub',xcen(i,:)'+rad);
                % Get precise partition of inCover.Theta using intersection
                Theta = inCover.Theta & boxTheta;
                % Instantiate the coverObj for each partition
                outCovers{i} = coverObj(Theta,xcen(i,:)',[],inCover.MPCi,boxTheta,inCover.T,rad);
            end
        end
    end % end methods
end % end classdef