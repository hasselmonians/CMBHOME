function state = vec2State(vec)
% Takes in a vector (nX1 or 1Xn) and reformats to a string, so can be used
% in MySQL statements.

vec = vec(:);

state = '(';

for i = 1:length(vec)
    state = [state vec(i) ',']; 
end

state = state(1:end-1);

state = [state ')'];


end