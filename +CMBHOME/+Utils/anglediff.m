function k = anglediff(varargin)
    p = inputParser;
    p.addRequired('A',@isnumeric);
    p.addOptional('B',NaN,@isnumeric);
    p.addParamValue('r',pi,@(x)isnumeric(x)&&((isscalar(x)&&x>0)||(numel(x)==2&&x(2)>x(1))));
    if nargin==3&&~any(cellfun(@(x)isequal(x,'r'),varargin))
        varargin = [varargin(1:2) {'r'} varargin(3)];
    end
    p.parse(varargin{:});
    
    r = p.Results.r;
    A = p.Results.A;
    B = p.Results.B;
    
    if isscalar(r)
       r = [-r r]; 
    end
    
    if ~all(isnan(B))
       k = A-B;
    else
       k = diff(A); 
    end
    
    k=(k-mean(r))/(range(r)/2);
    k(k<-2) = -mod(k(k<-2),2);
    k(k>2) = mod(k(k>2),2);
    k(k>1) = -2+k(k>1);
    k(k<=-1) = 2+k(k<=-1);
    k = k*(range(r)/2)+mean(r);
end