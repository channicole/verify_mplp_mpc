% Author:   Nicole Chan
% Created:  1/10/18
% Description: The 'tree' object is 
%
classdef tree < handle
    properties
        rootNode    % handle to root node
        size        % number of nodes in tree
        height      % height of tree
    end
    
    methods
        % Constructor
        function node = treeNode(root)
            if nargin ~= 1
                node.rootNode = treeNode();
            else
                node.rootNode = root;
            end
            node.size = 1;
            node.height = 0;
        end
        
        % Returns handle to root node
        function root = getRoot(node)
            root = node.rootNode;
        end
        
    end % end methods
end % end classdef