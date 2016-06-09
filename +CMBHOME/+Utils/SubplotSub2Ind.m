function ind = SubplotSub2Ind(size, pos)
% 2D subplot position to 1d index (for subplot)
%
% pos is the row, col (top left corner is 1,1)
% size is the [rows, cols] of subplot
%ind = CMBHOME.Utils.SubplotSub2Ind(size(mat),[row col])

M = size(2);

if any(size-pos < 0), error('POS cannot be larger than SIZE'); end

ind = M * (pos(1)-1) + pos(2);