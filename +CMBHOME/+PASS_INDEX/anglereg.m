function [ s,b ] = anglereg( x, theta, tolerance, maxiter )
% ANGLEREG - Linear-circular regression
%   Performes a linear-circular regression (Kempter et al, 2012)
%
%   [ S,B ] = ANGLEREG(X,THETA)
%   [ S,B ] = ANGLEREG(X,THETA,TOLERANCE)
%   [ S,B ] = ANGLEREG(X,THETA,TOLERANCE,MAXITER)
%
%   ARGUMENTS
%   * X: A vector of linear values
%   * THETA: A vector, the same size as X, containing circular values (in
%   radians)
%   * TOLERANCE (Optional): Default 0. Minimum step size persued in the
%   regression
%   * MAXITER (Optional): Default Inf. Maximum number of iterations for the
%   regression
%
%   RETURNS
%   * S: The slope in cycles per unit. To convert to radians, multiply by
%   2*pi
%   * B: The phase offset in radians.
%
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com)
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

step = 500;
rng = [-50*pi 50*pi];
n = length(x);

nn = 0;

while nn<maxiter&&range(rng)>tolerance
    s = unique(linspace(rng(1),rng(2),step)');

    [~,k]=max(sqrt((sum(cos(repmat(theta,[length(s) 1])-s*x),2)/n).^2+(sum(sin(repmat(theta,[length(s) 1])-s*x),2)/n).^2));
    rng = s(k)+[-2 2]*range(rng)/step;

    nn=nn+1;
end

s = s(k);

b = atan2(sum(sin(theta-2*pi*s*x)),sum(cos(theta-2*pi*s*x)));
