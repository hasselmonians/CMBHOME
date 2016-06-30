function tf = subStrCmp(C,S)
% function tf = subStrCmp(C,S)
%
% Test if string S (S can be a cell array of strings) exists in C (C can be
% cell array of strings).  If both are cell arrays, results will have size
% N x M where N is the length of cell array C and M is the length of cell
% array S.

if ~iscell(C), C = {C}; end

if ~iscell(S)
  
  tf = cellfun(@(c) ~isempty(c),strfind(C,S),'unif', 1);
  
else
  
  for i = 1:size(S,1)
    for j = 1:size(S,2)
      if isnan(S{i,j})
        tf(:,i,j) = false(size(C));
      else
        tf(:,i,j) = cellfun(@(c) ~isempty(c),strfind(C,S{i,j}),'unif', 1);
      end
    end
  end
end
