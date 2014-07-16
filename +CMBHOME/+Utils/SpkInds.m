function inds = SpkInds(t, spk)

% this function replaces spike_location, which used to incorrectly align
% spike times with video stamps and spit out spk_x, spk_y and spk_t. Now,
% we just have spk_ind, which are indices in structin.x, etc

count = hist(spk, t);

inds=[];

for i = 0:max(count)
    
    inds = cat(2, inds, find(count>i));
    
end

inds = sort(inds);

inds = inds(:);