function [ X ] = anglewrap( X,r)

if ~exist('r','var')
    r = [-pi pi];
end

if isscalar(r)
    r = [-r r];
end

X = (X-r(1))/range(r);
X(X<=0)=1+rem(X(X<0),1);
X(X>1)=mod(X(X>1),1);
X = X*range(r)+r(1);

end

