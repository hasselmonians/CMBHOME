function [ind,tf] = OverThresholdDetect(A, thresh, min_sep, min_length)
% Searches for continuous epochs for which A is at least thresh, merges those
% separated by less than min_sep indices, and returns ind, an Nx2 array of
% indices for continuous epochs over thresh
% ind = ThresholdDetect(A, thresh, min_sep, min_length);
% [ind, epochs] = ThresholdDetect(A, thresh, min_sep, min_length);

A = A(:);

tmp_ind = find(A<thresh); % all A under thresh to be thrown out

if all(A<thresh)
    
    ind = []; % none are over thresh

elseif ~isempty(tmp_ind);

    ind(:,1) = [1;tmp_ind+1];

    ind(:,2) = [tmp_ind-1; length(A)];

    ind(ind(:,2)-ind(:,1) < 0, :) = []; % remove all contiguous crossings
    
    if isempty(ind), ind = [tmp_ind(1), tmp_ind(end)]; end 
    
    % merge those close enough as per max_sep
    
    ind = reshape(ind', 1, numel(ind)); % reshape into vector like: start, stop, start2, stop2, start3...

    tmpdiff = diff(ind);
    tmpdiff(1:2:end) = min_sep+5; %make all inter-epoch times greater than max_sep
    
    mergeinds = find(tmpdiff<min_sep);

    mergeinds = cat(2, mergeinds, mergeinds+1);

    ind(mergeinds) = [];    % remove all overlapping epochs

    ind = reshape(ind, 2, numel(ind)/2)';

    % remove epochs less than min_length

    ind(ind(:,2) - ind(:,1)<min_length-1,:) = [];
    
else
    
    ind = [1 length(A)];
    
end


tf = zeros(length(A),1); % build tf logical vector

for i = 1:size(ind,1)
    
    tf(ind(i,1):ind(i,2)) = 1;
    
end

tf = logical(tf);