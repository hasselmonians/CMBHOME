function [ s,b ] = anglereg( x, theta, tolerance, maxiter )


if ~exist('tolerance','var')
    tolerance=0;
end

if ~exist('maxiter','var');
    maxiter=inf;
end

theta = mod(theta,2*pi);

if size(x,1)>size(x,2)
    x = x';
end
if size(theta,1)>size(theta,2)
    theta = theta';
end

step = 100;
rng = [-50*pi 50*pi];
n = length(x);

nn = 0;

while nn<maxiter&&range(rng)>tolerance
    
    s = unique(linspace(rng(1),rng(2),step)');
    
    %         plot(s,sqrt((sum(cos(repmat(theta,[length(s) 1])-s*x),2)/n).^2+(sum(sin(repmat(theta,[length(s) 1])-s*x),2)/n).^2));
    %         pause(0.5);
    
    [~,k]=max(sqrt((sum(cos(repmat(theta,[length(s) 1])-s*x),2)/n).^2+(sum(sin(repmat(theta,[length(s) 1])-s*x),2)/n).^2));
    rng = s(k)+[-2 2]*range(rng)/step;
    
    nn=nn+1;
end

s = s(k);

%%
rng = [0 2*pi];

nn=0;

while nn<maxiter&&range(rng)>tolerance
    b = unique(linspace(rng(1),rng(2),step)');
    [~,k] = min(sum((mod(repmat(s*x,[length(b) 1])+repmat(b,[1 n]),2*pi)-mod(repmat(theta,[length(b) 1]),2*pi)).^2,2));
    rng = b(k)+[-2 2]*range(rng)/step;
end

b = b(k);