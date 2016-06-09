function epochs = MergeEpochs2(epochs, mergetouching)
% merge epochs of arbirary precision without creating binary vector
%
% epochs = CMBHOME.Utils.MergeEpochs2(epochs, mergetouching)
%
%
% returned sorted
%
% andrew 10/13/2011
% Updated: Ehren Newman, 19th November 2012. Sped-up.

OE = epochs;
if ~exist('mergetouching', 'var'), mergetouching = 1; end

if size(epochs,1)==1, return; end

if any(epochs(:,2)-epochs(:,1)<0), error('epochs must be zero or positive windows of time like [tstart tstop]'); end

[~, inds] = sort(epochs(:,1));

if mergetouching
    
    op = @ge;
    
else
    
    op = @gt;
    
end

epochs = epochs(inds,:); % sort epochs by start time

i = 1;

if 1
  comps = op(epochs(1:end-1,2), epochs(2:end,1));
  
  while any(comps)
    
    epochsTmp(:,1) = [epochs(1,1); epochs(find(~comps)+1)];
    epochsTmp(:,2) = [epochs(~comps); epochs(end,2)];
    epochs = epochsTmp; clear epochsTmp;
    comps = op(epochs(1:end-1,2), epochs(2:end,1));
    
  end
  
else
  
  while i < size(epochs,1)
    
    if eval(['epochs(i,2)' op{1} 'epochs(i+1,1)'])
      
      epochs(i, 2) = epochs(i+1, 2);
      
      epochs(i+1, :) = [];
      
    else
      
      i = i+1;
      
    end
    
  end
  
end