function [ rho,p,s,b ] = kempter_lincirc( x,theta,varargin )
% ANGLEREG - Linear-circular correlation
%   Performes a linear-circular correlation (Kempter et al, 2012)
%
%   [ RHO,P,S,B ] = kempter_lincirc(X,THETA)
%   [ RHO,P,S,B ] = kempter_lincirc(X,THETA,S,B)
%
%   ARGUMENTS
%   * X: A vector of linear values
%   * THETA: A vector, the same size as X, containing circular values (in
%   radians)
%   * S: The presumed best slope in cycles/unit. If not entered, uses
%   anglereg
%   * S: The presumed best phase shift in radians. If not entered, uses
%   anglereg
%
%   RETURNS
%   * rho: The linear-circular correlation coefficient
%   * p: The statical significance of the correlation
%   * S: The slope in cycles per unit. To convert to radians, multiply by
%   2*pi
%   * B: The phase offset in radians.
%
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com)
import CMBHOME.PASS_INDEX.*;

if nargin==1, s=varargin{1}; varargin = {}; end;
if strcmp(class(x),'single')
    x = double(x);
end

goods = ~isnan(x)&~isnan(theta);
x = x(goods);
theta = theta(goods);
if ~exist('s','var')
    [s,b]=anglereg(x,theta,varargin{:});
end

n = length(x);

phi = mod(s*x,2*pi);
theta = mod(theta,2*pi);

phi_ = angle(sum(exp(1i*phi))/n);
theta_ = angle(sum(exp(1i*theta))/n);

rho = abs(sum(sin(theta-theta_).*sin(phi-phi_))/...
    sqrt(sum(sin(theta-theta_).^2)*sum(sin(phi-phi_).^2)))*sign(s);

lambda = @(i,j)n^-1*sum((sin(phi-phi_).^i).*(sin(theta-theta_).^j));
z = rho*sqrt(n*lambda(2,0)*lambda(0,2)/lambda(2,2));
p = 1-erf(abs(z)/sqrt(2));

end

