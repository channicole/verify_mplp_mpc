% % Copied directly from MATLAB's website with very minor change
% % UPDATED 11/30/17: now allows for single previous pointer and multiple
% % next pointers.
classdef dlnode < handle
   % dlnode A class to represent a doubly-linked node.
   % Link multiple dlnode objects together to create linked lists.
   properties
%       Data
   end
   properties(SetAccess = private)
      Next = dlnode.empty(0,0);
      Prev = dlnode.empty;
      parInd = [];
   end
   
   methods
      function node = dlnode()
%          % Construct a dlnode object
%          if nargin > 0
%             node.Data = Data;
%          end
      end
      
      function insertAfter(newNode, nodeBefore)
         % Insert newNode after nodeBefore.
         removeNode(newNode);
%          newNode.Next = nodeBefore.Next;
         newNode.Prev = nodeBefore;
%          if ~isempty(nodeBefore.Next)
%             nodeBefore.Next.Prev = newNode;
%          end
         nodeBefore.Next{end+1} = newNode;
         newNode.parInd = length(nodeBefore.Next);
      end
      
%       function insertBefore(newNode, nodeAfter)
%          % Insert newNode before nodeAfter.
%          removeNode(newNode);
%          newNode.Next = nodeAfter;
%          newNode.Prev = nodeAfter.Prev;
%          if ~isempty(nodeAfter.Prev)
%             nodeAfter.Prev.Next = newNode;
%          end
%          nodeAfter.Prev = newNode;
%       end
      
      % TODO: CHECK CORRECT
      function removeNode(node)
         % Remove a node from a linked list.
%          if ~isscalar(node)
%             error('Input must be scalar')
%          end
         prevNode = node.Prev;
         nextNode = node.Next;
         if ~isempty(prevNode) && ~isempty(nextNode)
            prevNode.Next{node.parInd} = nextNode;
         end
         if ~isempty(nextNode)
             for i=1:length(nextNode)
                next = nextNode{i};
                next.Prev = prevNode;
             end
         end
         node.Next = dlnode.empty;
         node.Prev = dlnode.empty;
      end
      
      % TODO: CHECK CORRECT
      function clearList(node)
         % Clear the list before
         % clearing list variable
         prev = node.Prev;
         next = node.Next;
         removeNode(node)
         while ~isempty(next) % FIX THIS TO WHILE LOOP
            node = next;
            next = node.Next;
            removeNode(node);
         end
         while ~isempty(prev)
            node = prev;
            prev = node.Prev;
            removeNode(node)
         end
      end
   end
   
   methods (Access = private)
      function delete(node)
         clearList(node)
      end
   end
end