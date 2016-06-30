function [fnameArray] = path2array(fname)
% function [fnameArray] = path2array(fname)
% 
% takes in a string representing the path to a file or folder and breaks it
% into a cell array in which each cell holds the portions of the filename
% that were between file seperator symbols ('/' or '\'). 
%
% This will help use a single set of paths that can be used on both windows
% and linux/mac based systems.
% 
% INPUTS ARGS:
%  fname - a string with the filename with path.
%
% OUTPUT ARGS:
%  fnameArray - a cell array with as many cells as there were directories
%             plus one for the name of the file
%

% 100519 eln

if iscell(fname)
	fnameArray = fname;
	return
end

fs = filesep;

fileBrks = strfind(fname,fs);

if isempty(fileBrks)
  if fs=='\'
    fs = '/';
  else
    fs = '\';
  end
  fileBrks = strfind(fname,fs);
  
  if isempty(fileBrks)
    fnameArray = {fname};
    return
  end
end

% remove fs symbols at the beginning of fname if it exists.
if any(fileBrks==1)
  fname(1) = [];
  fileBrks = strfind(fname,fs);
end

% add fs to end of fname if not as end marker
if fileBrks(end)~=length(fname)
  fname(end+1) = fs;
  fileBrks = strfind(fname,fs);
end

% chunk path into cell array
lastInd = 1;
for i = 1:length(fileBrks)
  fnameArray{i} = fname(lastInd:fileBrks(i)-1);
  lastInd = fileBrks(i)+1;
end
  