function [vals inds] = minMag(x,dim)
% utility to find value closest to zero
% function vals = minMag(x,[dim])

% eln 101104

if isempty(x)
	vals = []; inds = [];
	return;
end

  if exist('dim','var')
    [~,inds] = min(abs(x),[],dim);
    if isvector(inds)
      siz = size(x); siz2 = siz; siz2(dim) = [];
      inds2 = sub2ind(siz,1:max(siz2),inds(:)');
      vals = x(inds2);
    else
      for d = 1:ndims(x)
        if d==dim
          subs{d} = inds(:);
        else
          repSize = size(x); repSize(dim) = 1; repSize(d) = 1;
          subs{d} = reshape(repmat(shiftdim([1:size(x,d)]',1-d),repSize),[],1);
        end      
      end
      inds2 = sub2ind(size(x),subs{:});
      vals = x(inds2);
      outSize = size(x); outSize(dim)=1;
      vals = reshape(vals,outSize);
    end
  else
    % consider getting rid of this and just defining 'dim' when needed as
    % 1.
    [~,inds] = min(abs(x));
    if isvector(inds)
      inds2 = sub2ind(size(x),inds(:)',1:size(x,2));
      vals = x(inds2);
    end
  end
    
end

