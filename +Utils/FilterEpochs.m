function epoch = FilterEpochs(epoch, mindur, minsep)
%
% epoch = CMBHOME.Utils.Filterepoch(epoch, mindur, minsep)
%
% For an Nx2 array of epoch values (time or otherwise), merges those closer
% than 'minsep' and removes those less than 'mindur'
%
% andrew 6/29/2011

[~, inds] = sort(epoch(:,1)); % sort ascending

epoch = epoch(inds,:)'; % 2xN array 

seps = epoch(1,2:end) - epoch(2, 1:end-1); % check for minimum separation

seps = find(seps<minsep);

epoch = epoch(:);

epoch([seps*2, seps*2+1]) = [];

epoch = reshape(epoch, 2, length(epoch)/2)'; % back to Nx2 epoch format

durs = epoch(:,2) - epoch(:,1);

epoch(durs<mindur,:) = []; % remove epoch that are too short