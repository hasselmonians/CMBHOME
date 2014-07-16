function [ b ] = isempty( self )
% Returns boolean of whether the session has cells

b = builtin('isempty',self);

if ~b
    b=isempty(self.cells);
end


end

