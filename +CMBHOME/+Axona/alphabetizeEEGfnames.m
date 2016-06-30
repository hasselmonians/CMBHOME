function [fnamesSorted,inds] = alphabetizeEEGfnames(fnames)

if ~iscell(fnames)
  fnameLengths = arrayfun(@(x) length(x.name), fnames, 'unif', 1);
else
  fnameLengths = cellfun(@(x) length(x), fnames, 'unif', 1);
end
  lengths = unique(fnameLengths);

fnamesSorted = [];
inds = [];
for i = 1:length(lengths)
	fnamesSorted = [fnamesSorted; fnames(fnameLengths==lengths(i))];
	inds = [inds; find(fnameLengths==lengths(i))];
end
