function [ a ] = cell2mat2( b )
 % cell2mat, but rotates vectors where 90* off
if iscell(b)
    
    if ~(all(cellfun(@(x)size(x,1),b)==size(b,1))||all(cellfun(@(x)size(x,2),b)==size(b,2)))
        b = b';
    end      
    
    a=cell2mat(b);
else
    a=b;
end

end

