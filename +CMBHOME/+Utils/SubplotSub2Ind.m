function ind = SubplotSub2Ind(size, pos)
% ind = CMBHOME.Utils.SubplotSub2Ind(size(mat),[row col])
%
% pos is the row, col (top left corner is 1,1)
% size is the [rows, cols] of subplot
% andrew oct 13 2010

M = size(2);

if any(size-pos < 0), error('POS cannot be larger than SIZE'); end

ind = M * (pos(1)-1) + pos(2);