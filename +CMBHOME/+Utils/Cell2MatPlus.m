function matout = Cell2MatPlus(cellin)

% returns a matrix from cell cellin, where differently sized vectors are
% padded with NaNs

% andrew april 28 2010

if ~iscell(cellin)
    matout = cellin;
    return
end

max_el = max(cellfun(@numel, cellin));

cellout = cellfun(@(c) [c(:)', nan(1, max_el-numel(c))], cellin, 'UniformOutput', false);

matout = vertcat(cellout{:})';