function ov = struct2mat(si)
%
% Creates an output matrix for each field of the input structure. Assumes
% that the field is same type between all cells.

fnames = fieldnames(si);

[nr nc] = size(si);

for i = 1:length(fnames)
    
    tv = eval(['si(i).' fnames{i} ';']);
    
    if (~iscell(tv) & numel(tv) ==1 & ~isstr(tv)) %if can be a single double value, then do so
    
        eval(['ov.' fnames{i} '=NaN(' num2str(nr) ',' num2str(nc) ');'])
        for r = 1:nr
            r = num2str(r);
            for c = 1:nc
                c = num2str(c);
                eval(['ov.' fnames{i} '(' r ',' c ')=si(' r ',' c ').' fnames{i} ';']);
            end
        end   
    else % otherwise, make it a cell.
        eval(['ov.' fnames{i} '=cell(' num2str(nr) ',' num2str(nc) ');'])
        for r = 1:nr
            r = num2str(r);
            for c = 1:nc
                c = num2str(c);
                eval(['ov.' fnames{i} '{' r ',' c '}=si(' r ',' c ').' fnames{i} ';']);
            end
        end  
        
    end
end
            
