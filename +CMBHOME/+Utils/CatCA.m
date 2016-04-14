function matout = CatCA(cellin)
% matout = CMBHOME.Utils.CatCA(cellin)
%
% Converts a cell array of vectors 'cellin' of size MxN to a matrix padded by NaNs of
% size max_vector x N x M. This is useful for data structures in CMBHOME.
% If cellin is a matrix, matout = cellin. The reason M (epochs) gets pushed
% back to 3rd dim is that it is more often the case there are multiple
% cells, so the vectorization is for the cells.
%
% andrew 26 may 2010

if isnumeric(cellin), matout = cellin; return; end

%areallvectors = cellfun(@isvector, cellin); % check to see if they are all vectors

% if any(~areallvectors)
%     error('cellin must contain only vectors');
% end

mc=max(max(cellfun(@numel,cellin)));

cellin = cellfun(@(x) [x;nan(mc-numel(x),1)], cellin,'uni', false); % pad all with zeros

cellin = reshape(cellin, 1, size(cellin,2), size(cellin, 1)); % reshape cellin

matout = cell2mat(cellin);