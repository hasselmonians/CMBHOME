function [x,phi,a,phi_knot,rho_c,p] = corrC2L(x,phi) 
% Find the linear correlation between a linear and circular variable as in 
% R. Kempter 2012, Journal of Neuroscience Methods. 
%
% Inputs: x - the linear variable
%         phi - circular variable (must be in radians)
%           
% Output: x - Same as the input linear
%         phi - Same as input circular
%         a - The best-fit slope (radians per unit of x)
%         phi_knot - Best-fit intercept (radians)
%         rho_c - correlation coefficient
%         p - Significance value (unreliable for small datasets)
%
% [x,phi,slo,inter,rho,p] = utils.corrC2L(x,phi)

bads = isnan(x) | isnan(phi);
x(bads) = [];
phi(bads) = [];

n = length(phi);

f = @(a,phi,x) 1/sqrt(((1/n) * sum(cos(phi-2*pi*a*x))).^2 + ((1/n) * sum(sin(phi-2*pi*a*x))).^2);
[a R exitf] = fminsearch(@(a) f(a,phi,x),[0]);

phi_knot = atan2((sum(sin(phi - 2*pi*a*x))),(sum(cos(phi - 2*pi*a*x))));
ev = phi_knot+x*a; %ev = expected values from the fit.

%% rho is the correlation coefficient:
% Convert x into a circular variable:
theta = mod(2*pi*abs(a)*x,2*pi);
% Find Averages
theta_bar = atan2(sum(sin(theta)) , sum(cos(theta)));
phi_bar = atan2(sum(sin(phi)) , sum(cos(phi)));

num = sum(sin(phi-phi_bar).*sin(theta-theta_bar));
den = sum((sin(phi-phi_bar)).^2 ) .* sum((sin(theta-theta_bar)).^2);
rho_c = num/den;

l20 = lambda(2,0,x,phi,phi_bar,theta,theta_bar);
l02 = lambda(0,2,x,phi,phi_bar,theta,theta_bar);
l22 = lambda(2,2,x,phi,phi_bar,theta,theta_bar);

z = rho_c * sqrt(n*l20*l02/l22);
p = 1- erf(abs(z) / sqrt(2));

end


function lv = lambda(i,j,x,phi,phi_bar,theta,theta_bar)

lv = (1/length(phi))*sum((sin(phi-phi_bar)) .^i .* (sin(theta - theta_bar)) .^j);

end