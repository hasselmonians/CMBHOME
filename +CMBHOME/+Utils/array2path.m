function [fname] = array2path(fnameArray,replaceString)
% Recombines the strings stored in fnameArray 
% with the system appropriate
% filesep symbol to generate a valid path.
%
% INPUT ARGS:
%  fnameArray - A cell array in which each cell holds the names of folders
%  and the last cell may hold a file name.  
% 
% OUTPUT ARGS:
%  fname - a string with the path using the system appropriate file
%  seperators (\ or /).  The final character is never a file seperator.
%
% [fname] = array2path(fnameArray)

% eln 100519
cv=0;
fnameArray = fnameArray(~cellfun(@isempty,fnameArray));

fname = [];

if strcmp(fnameArray{1},'dropboxPath')
    if ~exist('replaceString','var'), replaceString = [];end
    if isempty(replaceString), replaceString = dropboxPath; end
    fnameArray{1} = replaceString; 
end

for i = 1:length(fnameArray)
  fname = [fname filesep fnameArray{i}];
end

if strcmp(fname(2),'~')
  fname = fname(2:end);
end

if fname(3)==':'
    fname = fname(2:end);
end


fname = strrep(fname,'//','/');
