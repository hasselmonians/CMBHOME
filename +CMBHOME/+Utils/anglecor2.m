function [ rho,p,s,b ] = anglecor2( x,theta,varargin )

import CMBHOME.Utils.*;

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

% if strcmp(class(rho),'single')
%     keyboard
% end

t = sign(rho) .* Inf;

k = (abs(rho) < 1);
t(k) = rho(k).*sqrt((n-2)./(1-rho(k).^2));

p = 2*tcdf(-abs(t),n-2);


end

%%

