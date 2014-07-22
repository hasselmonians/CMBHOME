function [f_cell, vel_dim] = VelocityRate(self, cel, vel_dim)
%[F, V] = root.VelocityRate(cel, vel_dim)
%
% Calculates the rate of firing for cell cel ([tetrode index, cell index]), 
% for bins of running speed.
%
% If vel_dim is specified, bin edges of running speed are vel_dim
% If vel_dim is not specified, bin edges of running speed are [0 max_vel]
% in increments of 1 cm sec^-1.
%
% Remember that root.spatial_scale sets cm/pixel.
%
% If multiple epochs are set in root.epoch, then all data from each epoch
% are concatinated and analysed as continuous
%
% root.epoch needs to be 1x2 vector
%
% andrew 17 may 2010

if ~exist('vel_dim', 'var')
    vel_dim = 0:self.spatial_scale^-1:max(self.vel);
end

if isempty(self.b_vel), disp('Warning, no user specified velocity assigned. The simple derivative thus used may not be optimal for this analysis'); end

import CMBHOME.Utils.*

self = MergeEpochs(self);

[vel, spk_vel] = ContinuizeEpochs(self.vel, self.spk_vel(cel)); % get running speed at spikes

tao_v = histc(vel, vel_dim) .* self.fs_video^-1; % time spent in each bin of running speed

spikes_v = histc(spk_vel, vel_dim);% spikes in each bin of running speed

f_cell = spikes_v ./ tao_v;

vel_dim = vel_dim(:);