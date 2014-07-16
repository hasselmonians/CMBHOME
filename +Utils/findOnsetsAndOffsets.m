function [OnOffs] = findOnsetsAndOffsets(boolVec)
% function [OnOffs] = findOnsetsAndOffsets(boolVec)
%
% Returns list of aligned start and stops of chunks of 1's in a vector of
% 1's and 0's.
%
% INPUTS
%  boolVec - vector of 1's and 0's 
%
% OUTPUTS
%  startEnds - Nx2 list of indices of the first and last 1's for the N
%              contiguous blocks of 1's.
%

% created 6/16/11 eln
% update 7/5/11 arb line 1
% update 7/6/11 eln line 26

boolVec = boolVec(:)';

starts = find(diff(boolVec)==1);
ends = find(diff(boolVec)==-1);

% if starts out going correct speed, add 1 to starts
if boolVec(1)
  starts = [0 starts];
end

% if finishes going correct speed, add final value to ends
if boolVec(end)
  ends = [ends length(boolVec)];
end

OnOffs = [starts(:)+1 ends(:)];