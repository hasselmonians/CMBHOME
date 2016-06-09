function inds = SpkInds(t, spk)
% Aligns two different timebases by using indices

count = hist(spk, t);

inds=[];

for i = 0:max(count)
    
    inds = cat(2, inds, find(count>i));
    
end

inds = sort(inds);

inds = inds(:);