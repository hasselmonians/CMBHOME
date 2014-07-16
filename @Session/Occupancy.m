function [occupancy, xdim, ydim] = Occupancy(self, xdim, ydim, mergeepochs, dim)
%
% Returns occupancy matrix of time spent in bins of dimension 3x3cm^2 be
% default. occupancy = root.Occupancy(xdim, ydim) where will return a matrix of
% dimensions xdim by ydim(:,2)
%
% Remember that root.spatial_scale indicates cm/pixel in your experiment.
%
% Occupancy is only returned for data within root.epoch. If multiple epochs
% are selected, they are all concatinated and one matrix is returned.
% [occupancy, xdim, ydim] = root.Occupancy returns occupancy matrix and x
% and y dimensions
%
% andrew bogaard 17 may 2010

import CMBHOME.Utils.*

if ~exist('mergeepochs', 'var'), mergeepochs = 1; end

if ~exist('dim', 'var'), dim = 3; end % cm width of each bin ( dim^2 cm^2 bins)

if mergeepochs
    self = MergeEpochs(self); 
    [x, y] = ContinuizeEpochs(self.x, self.y);
else
    x = self.x;
    y = self.y;
end

if ( ~exist('xdim', 'var') || ~exist('ydim', 'var') ) || (isempty(xdim) || isempty(ydim))
    if iscell(x)
        xdim = min(vertcat(x{:})):self.spatial_scale^-1*dim:max(vertcat(x{:})); %edges of x and y dimensions

        ydim = min(vertcat(y{:})):self.spatial_scale^-1*dim:max(vertcat(y{:}));
    else
        xdim = min(x):self.spatial_scale^-1*dim:max(x); %edges of x and y dimensions

        ydim = min(y):self.spatial_scale^-1*dim:max(y);
    end
    
end

if ~iscell(x)

    occupancy = hist3([x, y], 'Edges', {xdim, ydim}) / self.fs_video;

else

    occupancy = zeros(length(xdim), length(ydim), length(x));
    
    for i = 1:length(x)
        
        occupancy(:, :, i) = hist3([x{i}, y{i}], 'Edges', {xdim, ydim}) / self.fs_video;
        
    end
    
end