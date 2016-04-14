function rm = decodeRatemaps(a)

nnn=0;
for i = 1:length(a.ratemap)
    if ~isempty(a.ratemap{i})
        nnn = nnn+1;
    end
end

rm = cell(nnn,1); % Up till this point is so that we can preallocate



cv = 1;
for i =1:length(a.ratemap)
    if ~isempty(a.ratemap{i})
        rm{cv}=eval(a.ratemap{i});
        cv = cv+1;
    end
end

end
