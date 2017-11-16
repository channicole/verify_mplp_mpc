% Author:   Nicole Chan
% Created:  11/10/17
% Description: Checks the intersection between the reachtube (a 
% reachtubeObj object) and unsafe sets (stored in a function that unsafeSet
% points to). Returns true if the intersection is empty, and false 
% otherwise.
%
function safeflag=checkSafety(reachtube,unsafeSet)
    safeflag = 1;
    % Check each reachset in the reachtube
    for i = 1:length(reachtube)     % TODO: replace length() or add to class
        % Check each reachset against each unsafe set
        for j = 1:length(unsafeSet) % TODO: replace length()
            if ~isEmptySet(reachtube(i) & unsafeSet(j))
                safeflag = 0;
                return;
            end
        end
    end  
end