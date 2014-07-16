function k = anglemean(theta,r,weight)
    if exist('r','var')
       if isscalar(r), r=[-r r]; end;
       theta = 2*pi*(theta-r(1))/range(r)-pi;
    end
    
    if ~exist('weight','var')
        weight = ones(size(theta));
    end
    
    k = angle(sum(weight.*exp(1i*theta))/sum(weight));
    
    if exist('r','var')
       k = (k+pi)/(2*pi)*range(r)+r(1);
    end
end