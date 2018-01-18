% Author:   Nicole Chan
% Created:  1/10/18
% Description: The 'tree' object is 
%
classdef tree < handle
    properties
        rootNode@treeNode   % handle to root node
        size                % number of nodes in tree
        height              % height of tree
        currNode@treeNode   % handle to some node
    end
    
    methods
        % Constructor
        function tree = tree(root)
            if nargin ~= 1
                tree.rootNode = treeNode();
                tree.size = 0;
                tree.currNode = [];
            else
                tree.rootNode = treeNode(root);
                tree.size = 1;
                tree.currNode = tree.rootNode;
            end
            tree.height = 0;
        end
        
        % Returns handle to root node
        function root = getRoot(treeObj)
            root = treeObj.rootNode;
        end
        
%         % Updates the currNode handle
%         % TODO: enumerate direction type, determine what directions are
%         % necessary, determine how best to advance to sibling if used.
%         function [treeObj,treeObj.currNode] = stepCurrNode(treeObj,direction)
%             switch direction
%                 case 'parent'
%                     treeObj.currNode = treeObj.currNode.parent;
%                 case 'leftSibling'
%                     i = derp;
%                     treeObj.currNode = treeObj.currNode.parent.children{i};
%                 case 'rightSibling'
%                     treeObj.currNode = [];
%                 case 'leftChild'
%                     treeObj.currNode = treeObj.currNode.children{1};
%                 case 'rightChild'
%                     treeObj.currNode = treeObj.currNode.children{end};
%                 otherwise
%                     disp('Invalid input argument. currNode not updated.');
%             end  
%         end
    end % end methods
end % end classdef