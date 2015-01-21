function [angle_occupancy, theta] = DirectionalOccupancy(self, binsize, Continuize)
% O = root.DirectionalOccupancy;
% [O, theta] = root.DirectionalOccupancy(binsize, Continuize);
%
% Returns O, the directional occupancy of animal in Session root. If
% root.epoch is one epoch, O is a column vector. If root.epoch is more than
% one epoch, O is a matrix of column vectors of size M x N where M is the
% length of the number of bins as per binsize. N is size(root.epoch,1)
%
% Arguments:
% binsize -> positive number indicating binwidth in degrees
% Continuize -> 0 or 1 indicating whether to concatinate together multiple epochs
%
% Returns:
% O -> time in seconds in each bin
% theta -> bins in degrees
%
% andrew 28 may 2010

import CMBHOME.Utils.*

if ~exist('binsize', 'var'), binsize = 6; end

if ~exist('Continuize', 'var'), Continuize = 0; end

theta=-180+binsize/2:binsize:180-binsize/2;

theta = theta(:);

if Continuize
    
    self = MergeEpochs(self); % make sure epochs dont overlap so we dont double count
    
    headdir = ContinuizeEpochs(self.headdir); % make matrix of continuized epochs, multiple cells  (samples x cells)  

else
    
    headdir = self.headdir; % cell array of head directions
    
    headdir = CatCA(headdir); % make matrix (samples x cells x epochs)
    
end

angle_occupancy = zeros(length(theta), size(headdir, 2), size(headdir, 3));

for i = 1:size(headdir,3) % loop through epochs, this is as close to entirely vectorized as possible

    angle_occupancy(:,:,i) = hist(headdir(:,:,i), theta) / self.fs_video ;
    
end