function [R, mean_theta, S] = MeanResultant(self, cel)
% [R, mean_theta, S] = root.MeanResultant(cel); 
%
% Returns length, angle and circular variance pertaining to the mean
% resultant of directionality of neuron spiking.
%
% RETURNS:
% R, mean_theta and S. M x N array where N is the number of cells, and
% M is the number of epochs
%
% see equation 10, DIRECTIONAL STATISTICS. Gary L. Gaile & James E.
% Burt. 1980
%
% andrew 10 june 2010

[F, theta] = self.DirectionalTuningFcn(cel);

theta = deg2rad(theta(:));

theta = repmat(theta, [1, size(F,2), size(F,3)]);

mean_theta = atan2(sum( F .* sin(theta) , 1), sum( F .* cos(theta), 1)); % resultant angle

R = (cos(mean_theta) .* sum( F .* cos(theta), 1) + sin(mean_theta) .* sum( F .* sin(theta), 1)) ./ sum(F, 1); % mean resultant length

S = 1-R; % circular variance

R = shiftdim(R,1)';
S = shiftdim(S,1)';
mean_theta = shiftdim(mean_theta, 1)';