function self = FixDir(self, max_allowed_flips)
% Corrects errors in headdirection tracking
%
% Assumed root.b_headdir is in degrees
%
% Max allowed flips is the maximum number of consecutive lost samples
% (default = 5)
%
% Neuralynx head direction recordings often have missed samples, registered
% as either 0 or NaN. This script finds all missed samples and accidental
% flips of 180 degree (or any large flips > 162 degress) that persist for more than int
% 'max_allowed_flips' (5 samples) and interpolates the directional data
% in-between. Then we smooth the vector by convolving it with a gaussian of
% standard deviation of 3 bins, (or .1 sec at 30hz sampling freq).
%
% For all numerical samples (0 excluded) we continuize the data so
% that there are no more jumps at 360 and 0, and linearly interpolate all
% areas of bad data (either large angular velocities or 450s, NaNs, or
% zeros) that are 'max_allowed_flips' or less. Then smooth, and bring back
% to [-180,180).
%
% root = root.FixDir(max_allowed_flips);

if self.raw_headdir~=1
    %disp('It appears the head direction tracking data has already been affected. Reset root.raw_headdir=1 to override.');
    return;
end

    if isempty(self.b_headdir), disp('FixDir() Error: No values in b_headdir'); return; end

    thetain = self.b_headdir;
    
    bads = (thetain==0 | thetain==450);
    thetain(bads) = NaN;

    % circular derivative
    thetain=pi*thetain./180; thetaout = thetain; % convert thetain to radians

    delta_theta=diff(thetain(~isnan(thetain))); % take difference in discontinuous angles in radius

    delta_theta=atan2(sin(delta_theta), cos(delta_theta)); % derivative of circular data

    thetaout(~isnan(thetaout)) = [thetain(find(~isnan(thetain),1,'first')); cumsum(delta_theta) + thetain(find(~isnan(thetain),1,'first'))]; % continuized data

    thetaout = CMBHOME.Utils.RemoveJumpsAndSmooth(thetaout, 5, .5*pi);

    thetaout=atan2(sin(thetaout),cos(thetaout)); % bring back to radians bound by [-pi, pi)

    thetaout=(thetaout*180)./pi;

    self.b_headdir = thetaout;
    self.raw_headdir =  self.raw_headdir-1;

end