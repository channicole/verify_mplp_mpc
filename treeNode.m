% Author:   Nicole Chan
% Created:  1/10/18
% Description: The 'treeNode' object is 
%
classdef treeNode < handle
    properties
        data        % main node info
        parent      % handle to parent node
        children    % array of children node handles (for now, using custom implementation of stack)
    end
    
    methods
        % Constructor
        function node = treeNode(data,parent,children)
            if nargin < 3 || ~isa(children,'stackObj')
                children = stackObj();
            end
            if nargin < 2
                parent = [];
            end
            if nargin < 1
                data = [];
            end
            node.data = data;
            node.parent = parent;
            node.children = children;
        end
        
        % Check if this node is the root
        function flag = isRoot(node)
            flag = 0;
            if isempty(node.parent)
                flag = 1;
            end
        end
        
        % Check if this node is a leaf
        function flag = isLeaf(node)
            flag = 0;
            if isempty(node.children)
                flag = 1;
            end
        end
        
        % Update data
        function node = setValue(node,newData)
            node.data = newData;
        end
        
        % Add child
        function node = insertChild(node,child)
            node.children.push(child);
        end
        
        % Remove child; NOTE: for now, just removes last added child
        function node = removeChild(node)
            node.children.pop();
        end
    end % end methods
end % end classdef