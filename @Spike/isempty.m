function [ t ] = isempty( self )

% Returns whether the spike structure is empty
        t = cellfun(@isempty,reshape({self.ts},size(self)));
end

