function [ X ] = nanInterp( X, method)
% x_interpolated = CMBHOME.Utils.nanInterp(x)
%
% Interpolates all Nan's in the input vector, according to the second input
% argument "method". Defaults to nearest neighbor interpolation. 

% Bill 2013.01.08

if ~exist('method','var'), method = 'nearest'; end;

if all(isnan(X))
    X=zeros(size(X));
    return
elseif isscalar(X)||isempty(X)
    return
else
    
    if isnan(X(1))
        X(1:find(~isnan(X), 1 )-1) = X(find(~isnan(X), 1 ));
    end
    
    if isnan(X(end))
        X(find(~isnan(X), 1, 'last' )+1:end) = X(find(~isnan(X), 1, 'last' ));
    end
    
    
    temp = (1:length(X))';
    
    X = interp1(temp(~isnan(X)),X(~isnan(X)),temp,method);
end

end

