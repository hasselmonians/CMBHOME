function epochs = IntersectEpochs(root, epochs)
% epochs = root.IntersectEpochs(epochs)
%
% Returns list of epochs that are contained in both epochs lists { epochs1, epochs2 }
%
% INPUTS
%  epochs - cell array of epoch lists
%
% OUTPUTS
%  epochs - the common epochs between the two epochs 
%
% arb sept 20 2011 (Ehren's is commented out below

import CMBHOME.Utils.*

ts = root.b_ts;

tf_ts = zeros(length(ts),2); % boolean epochs

root.epoch = epochs{1}; % see where first epochs are

tf_ts(ContinuizeEpochs(root.ind),1) = 1;

root.epoch = epochs{2}; % see where second epochs are

tf_ts(ContinuizeEpochs(root.ind),2) = 1;

tf_ts = tf_ts(:,1) & tf_ts(:,2);

iEpochs = CMBHOME.Utils.findOnsetsAndOffsets(tf_ts);

epochs = root.b_ts(iEpochs);


% function [epochs] = IntersectEpochs(root,epochs)
% % function [epochs] = IntersectEpochs(epochs)
% %
% % Return list of epochs that are contained in both epoch lists
% %
% % INPUTS
% %  epochs - cell array of epoch lists
% %
% % OUTPUTS
% %  epochs - the common epochs between the two epochs
% 
% if ~iscell(epochs)
%     error('Epochs must be a cell array.');
% end
% 
% maxt = max(max(vertcat(epochs{:})));
% 
% % pre allocate indices
% inds = zeros(length(epochs),find(root.b_ts==maxt));
% 
% % Process each epoch list
% for j = 1:length(epochs)
% 
%     root.epoch = epochs{j};
%     
%     % epochs to be merged usually each have multiple epochs stored in cells
%     if iscell(root.ind)
%       
%       iEpochs = cellfun(@(c) [c(1) c(end)], root.ind, 'unif', 0);
%       
%       iEpochs = vertcat(iEpochs{:});
%     
%     else % if an epoch to be merged only has one epoch, root.ind won't be in a cell
%       
%       iEpochs = [root.ind(1) root.ind(end)];
%       
%     end
%       
%     for i = 1:size(iEpochs,1)
%     
%         inds(j,iEpochs(i,1):iEpochs(i,2)) = 1;   
%     
%     end
% 
% end
% clear epochs    
% 
% % now find common ground
% commonInds = all(inds,1);
% 
% iEpochs = CMBHOME.Utils.findOnsetsAndOffsets(commonInds);
% epochs(:,1) = root.b_ts(iEpochs(:,1));
% epochs(:,2) = root.b_ts(iEpochs(:,2));