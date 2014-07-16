function root = MergeEpochs(root)
%
% root = CMBHOME.Utils.MergeEpochs(root);
%
% Returns root with epochs that do not overlap. Touching epochs are merged.
%
% Assumes that all epochs are on the timescale of root.b_ts (all epochs
% have samples within them. any epochs without samples within them are
% deleted)
%
% v1. andrew
% v2. ehren
    
if ~strcmp(class(root), 'CMBHOME.Session'), error('MergeEpochs has been updated to require new syntax.'); end

if size(root.epoch,1)<2, return; end

root.epoch(cellfun(@isempty, root.ts, 'unif', 1), :) = []; % delete epochs that are empty

iEpochs = cellfun(@(c) [c(1) c(end)], root.ind, 'unif', 0);
    
iEpochs = vertcat(iEpochs{:});
    
inds = zeros(1,max(iEpochs(:)));

for i = 1:size(iEpochs,1)
    
    inds(iEpochs(i,1):iEpochs(i,2)) = 1;   
    
end

iEpochs = CMBHOME.Utils.findOnsetsAndOffsets(inds);

iEpochs(iEpochs(:,2)-iEpochs(:,1)<=0,:) = [];

root.epoch = root.b_ts(iEpochs);