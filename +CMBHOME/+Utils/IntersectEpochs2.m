function epochs = IntersectEpochs2(varargin)
% epochs = IntersectEpochs2(epochs1, epochs2, epochs3...)
%
% Must pass at least
%
% Finds the intersection of all sets of epochs.
%
% andrew october 14 2011

import CMBHOME.Utils.*

if length(varargin)<2, error('you must pass at least 2 input arguments of type ''epoch'''); end

if any(cellfun(@isempty, varargin, 'unif', 1)), epochs = []; return; end % if any epochs are empty, there is no intersection

for i = 2:length(varargin)
   
    varargin{1} = MergeEpochs2(varargin{1}, 1);
    
		varargin{i} = MergeEpochs2(varargin{i}, 1);
		
    e_tmp = SelectEpochs(varargin{1}, varargin{i}, 1);
		
		e2_tmp = SelectEpochs(varargin{i}, varargin{1}, 1);
    
		varargin{1} = vertcat(e_tmp, e2_tmp);
		
end

epochs = varargin{1};