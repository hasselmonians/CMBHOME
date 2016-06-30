function [relPath] = relativePath(basePath,newPath)
% function [relPath] = relativePath(basePath,newPath)
%
% Returns the relative path from 'basePath' to 'newPath'

% eln 100415

% compute where all the file seperators are in the basepath
bp = fullfile(basePath,'.');
bpDirMrks = ismember(bp,filesep);
bpDirMrks = cumsum(bpDirMrks) .* bpDirMrks;

% there should be at least some directory structure
% throw an error if none are found
if max(bpDirMrks)==0
  error('No fileseparators found in basePath.\n');
end

% loop through directories in the basePath
for brk = 1:max(bpDirMrks)
  % build up basePath to check extent of overlap
  currBp = bp(1:find(bpDirMrks==brk));

  % quit now if base path is longer than newPath
  if length(currBp)<=length(newPath)

    % quit now if basePath and newPath no longer overlap
    if ~strcmp(newPath(1:length(currBp)),currBp)

      % set brk back to the last one that worked
      brk = brk - 1;
      break
    end
  else

    % set brk back to the last one that worked
    brk = brk - 1;
    break
  end
end % for brk

% back up directory structure as needed
relPath = repmat(['..',filesep],1,max(bpDirMrks)-brk);

% take only the part of the base dir that still overlapped
currBp = bp(1:find(bpDirMrks==brk));

% assemble the relative path 
relPath = [relPath,newPath(length(currBp)+1:end)];

