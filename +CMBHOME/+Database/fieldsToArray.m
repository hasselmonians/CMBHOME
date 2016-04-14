function ov = fieldsToArray(q)
% Takes in a 1x1 Struct with non-scaler fields and outputs nx1 struct with
% scaler (or cell) fields for each row. 

% wchapman 2013.03.06

fn = fieldnames(q);

lv = eval(['length(q.' fn{1} ')']);

for i = lv:-1:1
    for k = 1:length(fn)
        eval(['ov(i).' fn{k} '= q.' fn{k} '(i);']);
    end
end



end