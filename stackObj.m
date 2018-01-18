% Author:   Nicole Chan
% Created:  10/6/17
% Description: Implementation for a stack. Dynamically resizes with each
% push and pop.
%
classdef stackObj < handle
    properties
        items = stackObj.empty;    % objects in stack
        last  = [];             % index of last object in stack
    end
    
    methods
        function s = stackObj(obj)
            if nargin == 1
                if iscell(obj)
                    dim = numel(obj);
                    s.items = obj;
                    s.last = dim;
                else
                    s.items = {obj};
                    s.last = 1;
                end
            else
                s.items = {};
                s.last = 0;
            end
        end
        function push(s,obj)
            if iscell(obj)
              dim = numel(obj);
              s.items(s.last+1:s.last+dim) = obj;
              s.last = s.last+dim;
            else
              s.last = s.last+1;
              s.items(s.last) = {obj};
            end
        end
        function pop(s)
            if s.last~=0
                s.items(s.last) = []; % use 's.items{s.last} = []' if not deleting cell
                s.last = s.last-1;
            end
        end
        function y = isempty(s)
          y = (s.last == 0);
        end
        function y = size(s)
          y = s.last;
        end
        function y = getLast(s)
          if s.last~=0
              y = s.items{s.last}; 
          else
              y = [];
          end
        end
    end
end