% Author:   Nicole Chan
% Created:  1/11/18
% Description: Computes the reachtube of the given initial set of states
% over the given time horizon for a hybrid system.
%
% When this function is called, the input 'cover' should already correspond
% to the dynamics of no more than one discrete mode. (TODO: add check if
% desired).
%
function cover=computeReach(cover,in)
    if nargin~=2
        error('Not enough input arguments.');
    end
    %% Start building tree of reachObj nodes (these nodes are eventually concatenated into coverObj.Reach)
    RTree = tree(cover.data.getReachtube); % Note RTree.currNode is set to the root node
    
    while ~isempty(RTree.currNode)
        % Compute Post (only leaf nodes do NOT need Post)
        if RTree.currNode.data.needsPost() % RTree.currNode.data should be type 'reachtubeObj'
            tf = min(RTree.currNode.data.T(end),RTree.currNode.data.T(1)+in.vPar.deltaStep);
            [outReach,outSafe] = computePost(RTree.currNode.data.Theta,...
                RTree.currNode.data.T(1:find(RTree.currNode.data.T==tf,1)),...
                in.flowEq,in.MPCsol,in.modelJac,in.unsafeStates,in.vPar.epsilonConst,in.vPar.LipConst);
            
            % If found a possibly unsafe region, stop reachtube computation
            if ~outSafe
                cover.data.setSafeFlag(safetyEnum.NotSafe);
                % Update currNode
                RTree.currNode = treeNode.empty();
                
            % If safe, continue computing reachtube
            elseif outReach.T(end) < RTree.currNode.data.T(end)
                % Update node with computed sub-reachtube
                [~,~,newReach] = outReach.getProperties();
                RTree.currNode.data.updateReach(newReach);
                
                % Add nodes/sets of initial sets for next interval of time
                newT = RTree.currNode.data.T(find(...
                    RTree.currNode.data.T == outReach.T(end),1):end);
                [newx0,newrad,newBall,newTheta] = poly2ball(newReach(end),in.MPCsol.Pn);
                for i=1:size(newrad,1)
                    RTree.currNode.insertChild(treeNode(reachtubeObj(...
                        newTheta(i),newx0(i,:),[],[],newBall(i),newT,newrad(i,:))));
                end
                
                % Set currNode to point to the left-most child
                RTree.currNode = RTree.currNode.children.getItem(1); %RTree.currNode.children.items{1}
                
            % If safe and completed a branch of the reachtube computation
            else
                % If a left sibling exists
                lSib = RTree.currNode.parent.children.getItem(1);
                if lSib~=RTree.currNode
                    if lSib.data.needsPost()
                        error('Left sibling should have Post computed already. Derp.');
                    end
                    % Combine reachsets computed in these sibling nodes
                    RTree.currNode.data.reachUnion(lSib.data);
                    % Pop left sibling
                    RTree.currNode.parent.removeChild(1); % if lSib is referencing this node, then we might need to clear lSib or Matlab might not delete the object itself
                end
            end
            
        % If a leaf node and right sibling exists, update currNode to right sibling
        elseif RTree.currNode.parent.children.size > 1
            RTree.currNode = RTree.currNode.parent.children.getItem(2); % A leaf node at this point should always be the left-most sibling; hopefully this is a handle to the correct node and not a copy; Try "RTree.currNode.parent.children.items{2}"
        
        % If a leaf node and no siblings, update currNode to parent
        elseif ~isempty(RTree.currNode.parent)
            % Combine reachsets with parent
            RTree.currNode.data.reachUnion(RTree.currNode.parent.data);
            % Update currNode to parent
            RTree.currNode = RTree.currNode.parent;
            % Delete child
            RTree.currNode.removeChild(1);
            % % SANITY CHECK % %
            if RTree.currNode.children.size ~= 0
                error('Child not correctly removed or siblings not corectly handled.');
            end
            % % % % % % % % % % %
        % No parent means root has been reached and reachtube computation is complete
        else
            % % SANITY CHECK % %
            if RTree.currNode~=RTree.rootNode || RTree.size~=0
                error('Motherlover.');
            end
            % % % % % % % % % % %
            % Update currNode now that tree search is complete
            RTree.currNode = treeNode.empty();
        end
    end

end