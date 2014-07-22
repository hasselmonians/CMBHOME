function [A, B, phi, x, y] = EllipseFromCoef(cs, XY)
%
% [majoraxis, minoraxis, phi_radians, vecx, vecy] = EllipseFromCoef(cs, XY);
%
% XY are the points on the ellipse. this is a cheater because I dont feel
% like solving for an acceptable range for y(x) from the coefficients
%
% takes the coffecients from EllipseDirectFit and turns them into
% meaningful attributes of an ellipse

if length(cs)<5
  [A, B, phi, x, y] = deal(nan);
  return
end

x = -3*max(XY(:)):.1:3*max(XY(:)); % build vectors to plot the ellipse

y = (-cs(2).*x - cs(5) + sqrt((cs(2).*x+cs(5)).^2 - 4*cs(3)*(cs(6)+cs(4).*x+cs(1).*x.^2))) / (2*cs(3));
y2 = (-cs(2).*x - cs(5) - sqrt((cs(2).*x+cs(5)).^2 - 4*cs(3)*(cs(6)+cs(4).*x+cs(1).*x.^2))) / (2*cs(3));

x = [x(imag(y)==0), fliplr(x(imag(y2)==0))];
if isempty(x),   [A, B, phi, x, y] = deal(nan);  return, end
x = [x, x(1)];

y = [y(imag(y)==0), fliplr(y2(imag(y2)==0))];
if isempty(y),   [A, B, phi, x, y] = deal(nan);  return, end
y = [y, y(1)];

cs = cs'.*[1 .5 1 .5 .5 1];

sa = [sqrt(2*(cs(1)*cs(5)^2 + cs(3)*cs(4)^2 + cs(6)*cs(2)^2 - 2*cs(2)*cs(4)*cs(5) - cs(1)*cs(3)*cs(6)) / ...
    ((cs(2)^2-(cs(1)*cs(3)))*(sqrt((cs(1)-cs(3))^2+4*cs(2)^2) - (cs(1)+cs(3)))));... % semiax2

    sqrt(2*(cs(1)*cs(5)^2 +cs(3)*cs(4)^2 + cs(6)*cs(2)^2 - 2*cs(2)*cs(4)*cs(5) - cs(1)*cs(3)*cs(6)) / ...
    ((cs(2)^2-(cs(1)*cs(3)))*(-sqrt((cs(1)-cs(3))^2+4*cs(2)^2) - (cs(1)+cs(3)))))]; % semiax1   

if cs(1)>cs(3) & cs(2)==0, phi = .5*pi; % solve for angle of major axis
elseif cs(2)==0 & cs(1)<cs(3), phi = 0;
elseif cs(2)~=0 & cs(1)<cs(3), phi = .5 * acot((cs(1)-cs(3)) / (2*cs(2)));
elseif cs(2)~=0 & cs(1)>cs(3), phi = pi/2 + .5 * acot((cs(1)-cs(3)) / (2*cs(2)));
else phi = 0; end

if diff(sa)>0
    phi = pi/2 + phi;
    A = sa(2); % major axis
    B = sa(1); % minor axis
else
    A = sa(1);
    B = sa(2);
end

if A<B, error('an assumption about the ellipse equations was wrong.'); end

end

