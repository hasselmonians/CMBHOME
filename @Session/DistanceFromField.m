function [R, dR, field, rm, peaks] = DistanceFromField(self, cel, strict, plot, peaks)
% [R, dR, field, rate_map, field_centers] = root.DistanceFromField(cel)
%
% returns the minimum distance from field centers. corresponds to root.x
% and root.y. cannot do multiple epochs yet.
%
% may program this so dR is derived from the root.vel vector in case the
% kalman is there. for now, I will convolve with a gaussain of std 2bins to
% get rid of the zeros
%
% ARGUMENTS
%   cel             2 element vector, (1)=tetrode (2)=cell
%   strict          0 or (1). If 'strict' is set to 0, then we retain fields on the edge of the ratemap if the peak firing
%                   rate is at least 50% the max firing rate in the rate map
%   plot            (0) or 1. If plot is set, we plot the ratemap and the
%                   detected peaks
%   peaks           Nx2 array of where the peaks are in the ratemap
%                   [X(:),Y(:)]. If not initialized, the function will
%                   detect peaks.
%
% RETURNS
%   R               distance to nearest field (pixels)
%   dR              direction of movement relative to closest field (neg is
%                   inbound)
%   field           each field has a unique ID. this is that.

if ~exist('strict', 'var'), strict = 1; end

if ~exist('plot', 'var'), plot = 0; end

R = []; % initialize variables
dR = [];
field = [];

import CMBHOME.Utils.*

[rate_map, xdim, ydim, occ] = self.RateMap(cel, 'continuize_epochs', 1, 'supress_plot', 1, 'std_smooth_kernel', 9, 'binside', 4, 'omit_noocc', 0);

if ~exist('peaks', 'var')

    [f, inds] = extrema2(rate_map); % peaks in rate map

    inds = OmitEdgeFields(rate_map, inds, f, strict);

    [Y, X] = ind2sub(size(rate_map), inds); %inds of peaks, I=rows J=columns

    X = xdim(X);
    Y = ydim(Y);
else 
    X = peaks(:,1)';
    Y = peaks(:,2)';
end

if plot
    imagesc(xdim, ydim, rate_map), axis off, set(gca, 'Ydir', 'normal'), hold on, axis equal
    xs = xlim;
    ys = ylim;
        
    colormap(gray)
    scatter(X, Y, '.g', 'sizedata', 50);
end

if isempty(X), return; end % no fields

% keyboard
R = zeros(length(self.x), 1); % initialize distance from field
field = zeros(length(self.x), 1);

for i = 1:1000:length(self.x) % go in 1000 sample chunks
    
    i2 = min([999+i, length(self.x)]);
    
    % this Nx1xM matrix is the distance from each point (N) to the closest fields (M)
    tmp = sqrt(sum((repmat([self.x(i:i2), self.y(i:i2)], [1 1 length(X)]) - repmat([permute(X, [3 1 2]), permute(Y, [3 1 2])], i2-i+1, 1)).^2, 2));

    [R(i:i2), I] = min(tmp, [], 3);

    field(i:i2) = I; % field indices
    
end

dR = [0; diff(R)];

dR = conv(dR, pdf('normal', -5:5, 0, 2), 'same');

rm.ratemap = rate_map;
rm.xdim = xdim;
rm.ydim = ydim;

peaks = [X(:), Y(:)];

%%MAKE MOVIE%%%%%%%%%%%%%%%%%%%%%%%%%
%
% h = figure('Position', [100 100 1000 400])
% axis tight
% 
% % Record the movie
% f_ind = 1;
% 
% for j = 51:5:500
%  
%     subplot(1, 2, 1)
%     
%     scatter(X, Y), ylim([ydim(1) ydim(end)]), xlim([xdim(1) xdim(end)]), hold on
%     
%     plot(self.x(j-50:j), self.y(j-50:j)), hold off
%     
%     subplot(1, 2, 2), bar(R(j)), ylim([0 max(R)])
%     
%     F(f_ind) = getframe(h);
% 
%     f_ind = f_ind+1;
% end
% % Play the movie ten times
% movie(F,10)
% 
% keyboard

function inds = OmitEdgeFields(rm, inds, f, strict)
% look at all the fields in the rate map (the largest contiguous part of
% the ratemap), and checks that none of the detected peaks are along the
% edge. does this by checking that no more than 1 out of range pixels are
% touching
%
% if 'strict' is set to 0, then we retain edge fields if the peak firing
% rate is at least 50% the max firing rate in the rate map, unless is it a
% "peninsula'
%
% also takes away fields that are less than half the height of the highest
% peak.

[I, J] = ind2sub(size(rm), inds); % find peak positions

if strict==1 | strict==0 % if strict is -1, then just filter by peak firing rate

    todelete = [];
    for k = 1:length(I) % iterate through all peaks

        n_unocc = (3*(I(k)==1 | I(k)==size(rm,1))+3*(J(k)==1 | J(k)==size(rm,2))); % number of unoccupied spces near (yes, we double count edges in the case of corners, but thats OK, because corners should be removed)

        if ~strict & n_unocc<5 % if we are not being strict, remove only corners
            if f(k)>.5*max(f)
                continue;
            end
        end

        if n_unocc>1
            todelete = [todelete, k];
        end

    end

    inds(todelete) = [];
    f(todelete) = []; 

end

inds(f<mean(f)-std(f)) = []; % remove peaks that less than the mean - 1std
