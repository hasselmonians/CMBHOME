function [R P b] = plot_velocityrate(self, cel, params)
% root.plot_velocityrate(cel)
% root.plot_velocityrate(cel, params)
%
% Plots the firing frequency of cell cel against animal running speed for
% root.epoch (if multiple, they are all concatinated. See root.VelocityRate
% for algorithm).
%
% Bins of velocity from 0cm/sec:1cm/sec:35 cm/sec (as per
% root.spatial_scale pixels/cm)
%
% ARGUMENTS
%
% cel -> vector with two elements: cel(1) = tetrode index,  cel(2) = cell
% index
%
% params (optional) -> 0 or 1 element vector indicates whether or not to plot a
% linear regression with correlation coefficient 
%
% andrew 10 june 2010

import CMBHOME.Utils.*

R = [];
P = [];
b = [];

if ~exist('params', 'var'), params = 1; end
  
vel_dim = 0:self.spatial_scale^-1:self.spatial_scale^-1*35; % go from 0 to 35cm/s using root.spatial_scale 1 cm at a time

[F, V] = self.VelocityRate(cel, vel_dim);

plot(V*self.spatial_scale, F, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 15, 'MarkerEdgeColor', 'k'), xlabel('Running Speed (cm/s)'), ylabel('Firing Frequency'), title('Firing Frequency v. Running Speed');

[Fn, R, P] = LinearRegression(V(~isnan(F)), F(~isnan(F)));

b = diff(Fn) ./ (diff(V(~isnan(F)))*self.spatial_scale);

b = b(1);

if params(1)
    
    hold on
            
    line(V(~isnan(F))*self.spatial_scale, Fn, 'Color', 'r', 'LineWidth', 1.5), hold on
    
    ys = ylim;
    xs = xlim;
    
    text(xs(2) - .9*diff(xs), ys(2)-.1*diff(ys), ['R = ' num2str(R, 3)], 'FontSize', 8)
    text(xs(2) - .9*diff(xs), ys(2)-.2*diff(ys), ['slope = ' num2str(b, 3)], 'FontSize', 8)
    
end

end
