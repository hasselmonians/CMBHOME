function [gridness, phi_v_cor, auto_corr, d]  = Gridness2(self, cel, varargin)
% Calculates the original gridness score of 'cel'. Continuizes all epochs.
%
% ARGUMENTS
%
%   cel             1x2 vector indicating tetrode index and cell index
%
% RETURNS
%
%   gridness        a number indicating gridness score
%   periodicity     vector of the autocorrelogram correlation as a function of
%                   rotation (see ALGORITHM below)
%   auto_corr_rm    matrix of the rate map autocorrelogram with donut cut
%                   out of it. There is a second matrix along the third
%                   dimension indicating where cuts were made. This can be
%                   used in the AlphaData property of an image plot
%
% OPTIONAL PARAMETERS
%   
%   xdim                vector of bin edges along x dimension
%   ydim                vector of bin edges along y dimension
%   supress_plot        0 or 1 (1). If 0, plots autocorrelation, the cut
%                       donut, and the periodicity of the donut
%   std_smooth_kernel   STD of the gaussian kernel to smooth the rate map  
%   binside             The length in cm of the side of a bin when
%                       calculating the rate map
%   rotate_inc          (3 degrees). Change the angle increment for which
%                       we calculate the correlation between rotated ac's
%   autocorr            you can provide an autocorrelation to run gridness
%                       on like this: CMBHOME.Session.Gridness([], [], 'autocorr', AC); 
%
% DEPENDS ON
%
%   image processessing toolbox (edge and padarray)
%
%   CMBHOME.Utils (extrema2 and others)
%
% ALGORITHM
%
% 1.    calculates the rate map, smooths it, sets no occupancy to zero
% 2.    gets the autocorrelation of the rate map
% 3.    cuts out the center peak from the autocorrelation, and finds the six
%       surrounding peaks, if they all exist, and cuts them out to make a
%       donut
% 4.    rotate the donut 3 (rotate inc) degrees at a time to 360, calculating
%       correlation at each step
% 5.    gridness score is the difference between the minimum peak at either
%       60 or 120 degrees and the maximum trough at either 30, 90 or 150
%       degrees
%
% [gridness, periodicity, auto_corr_rm] = root.Gridness(cel)
% [gridness, periodicity, auto_corr_rm] = root.Gridness(cel, params)
p = inputParser;

p.addRequired('self')
p.addRequired('cel', @isnumeric)
p.addParamValue('xdim', [], @isnumeric);
p.addParamValue('ydim', [], @isnumeric);
p.addParamValue('clims', [], @(c) numel(c)==1);
p.addParamValue('continuize_epochs', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('supress_plot', 1, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('figure_handle', [], @(c) numel(c)==1);
p.addParamValue('std_smooth_kernel', 2, @isnumeric);
p.addParamValue('binside', 3, @isnumeric)
p.addParamValue('rotate_inc', 3, @isnumeric)
p.addParamValue('autocorr', [], @isnumeric)

p.parse(self, cel, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
xdim = p.Results.xdim;
ydim = p.Results.ydim;
clims = p.Results.clims;
continuize_epochs = p.Results.continuize_epochs;
supress_plot = p.Results.supress_plot;
figure_handle = p.Results.figure_handle;
std_smooth_kernel = p.Results.std_smooth_kernel;
binside = p.Results.binside;
rotate_inc = p.Results.rotate_inc;
autocorr = p.Results.autocorr;

import CMBHOME.Utils.* % we need the super cool extrema2 function

if ~isempty(self)
    if size(self.epoch, 1)>1 && continuize_epochs==0 % use some recursion to do multiple epochs

        [gridness, phi_v_cor, auto_corr] = MultipleEpochs(self, xdim, ydim, binside, std_smooth_kernel,...
                                            rotate_inc, cel);

        if ~supress_plot, PlotIt(auto_corr, phi_v_cor, gridness); end

        return

    end
end

gridness = NaN;
phi_v_cor = nan(min(180/rotate_inc)+1, 2);

rot_ind_mins = round([30/rotate_inc, 90/rotate_inc, 150/rotate_inc])+1; % indexes in rotated ac to check

rot_ind_maxs = round([60/rotate_inc, 120/rotate_inc])+1;

if isempty(autocorr)

    [rate_map, ~, ~, ~, rmmask] = self.RateMap(cel, 'supress_plot', 1, 'xdim', xdim, 'ydim', ydim, 'binside', binside, 'std_smooth_kernel', std_smooth_kernel, 'continuize_epochs', 1); % step 1

    autocorr = moserac(rate_map); % step 2: get autocorrelogram, make sure its a square image for later processing...
    
else
    
    rate_map = [];
    
    rmmask = [];
    
end

auto_corr = MakeSquare(autocorr); % make autocorrelation matrix square

[auto_corr, d, mask, no_center_peak, thresh] = RemoveCenterPeak(auto_corr);

if no_center_peak, disp('No center peak'), return, end

r = ceil(d)+1:floor(size(auto_corr, 1)/2); % initialize variables
cell_auto_corr = cell(length(r), 1);
cell_mask = cell(length(r), 1);
gridness = nan(length(r), 1);
phi_v_cor = nan(min(180/rotate_inc)+1, 2*length(r)); 

for i = 1:length(r)

    [cell_auto_corr{i}, cell_mask{i}] = RemoveOuterCircle(auto_corr, r(i), mask);

    [~, ~, phi_v_cor(:,[(i-1)*2+1 (i-1)*2+2])] = AutoCorrRotation(rot90(cell_auto_corr{i}, 3), cell_auto_corr{i}, 'cut_circle', 0, 'supress_plot', 1, 'rotate_inc', rotate_inc);

    gridness(i) = min(phi_v_cor(rot_ind_maxs, (i-1)*2+2))-max(phi_v_cor(rot_ind_mins, (i-1)*2+2));
    
end

[gridness, ind] = max(gridness); % take max grid score

auto_corr = cell_auto_corr{ind};

mask = cell_mask{ind};

phi_v_cor = phi_v_cor(:,[(ind-1)*2+1, (ind-1)*2+2]);

phi_v_cor(:,1) = phi_v_cor(:,1) + 90;

auto_corr(:,:,2) = ~mask;

if ~supress_plot, PlotIt(auto_corr, phi_v_cor, gridness); end

end

function PlotIt(auto_corr, phi_v_cor, gridness)

figure

n_epochs = length(gridness);

for i = 1:n_epochs
    
    if iscell(auto_corr) 
        subplot(n_epochs, 2, (i-1)*2+1), imagesc(auto_corr{i}(:,:,1)), axis square, axis off
    else
        subplot(n_epochs, 2, (i-1)*2+1), imagesc(auto_corr(:,:,1)), axis square, axis off
    end
    
    subplot(n_epochs, 2, (i-1)*2+2), plot(phi_v_cor(:,(i-1)*2+1), phi_v_cor(:,(i-1)*2+2));
    
end

end

function [auto_corr, mask] = RemoveOuterCircle(auto_corr, r, mask)
%
% auto_corr         matrix of autocorrelation
% d                 half the width of peaks (radius of peak edge)
% mask              marks where cuts in auto_corr were made
% thresh            threshold for peaks

s = size(auto_corr, 1)-1; % square matrix

[col, row] = meshgrid(-s/2:s/2, -s/2:s/2); % indices on grid

auto_corr(sqrt(col.^2 + row.^2) > r) = 0;

mask(sqrt(col.^2 + row.^2) > r) = 1;

d_s = round((s+1 - r*2)/2);

auto_corr = auto_corr(d_s:end-d_s, d_s:end-d_s); % reshape autocorrelation

mask = mask(d_s:end-d_s, d_s:end-d_s); % reshape mask

end

% function [matout, d, mask, no_center_peak, thresh] = RemoveCenterPeak(matin)
% % [centercut, 2stdmajoraxis] = RemoveCenterPeak(matin);
% %
% % Removes center peak of autocorrelogram by finding point furthest from
% % center that is half the peak (which is == 1 if using moserac).
% %
% % The above is the radius of center cutout. The since we are using an
% % autocorr, the peak will be either round or an oblong shape. The ratio of
% % the half width at the minor axis to the major axis could be the skewedness of
% % fields to be added in a future version.
% %
% % andrew nov 10 2010
% 
% s = size(matin, 1)-1; % square matrix
% 
% no_center_peak = 0;
% 
% thresh = .1; % use sargolini et al threshold for peaks
% 
% [col, row] = meshgrid(-s/2:s/2, -s/2:s/2);
% % 
% % peaks = sqrt(col(matin >= thresh).^2 + row(matin>= thresh).^2);
% % 
% % d = sort(peaks(:));
% % 
% % d = d(find(diff(d)>1, 1, 'first')); % farthest radius of center peak
% 
% % edge_inds = edge(matin>=.1);
% % 
% % d = min(sqrt(col(edge_inds).^2 + row(edge_inds).^2));
% 
% img_prop = regionprops(matin>thresh, {'centroid', 'minoraxislength', 'majoraxislength', 'FilledArea'});
% 
% [~, ind] = min(sum((vertcat(img_prop(:).Centroid)-repmat((s+1)/2, length(img_prop), 2)).^2, 2)); % ind to center region
% 
% d = img_prop(ind).MinorAxisLength/2;
% 
% %d =  img_prop(ind).FilledArea * (img_prop(ind).MinorAxisLength/2) / ((img_prop(ind).MinorAxisLength/2)^2*pi);
% 
% if isempty(d) || 2*d>=.8 * (s+1) % if center peak is not a discrete thing, redefine thresh and look for closest non-thresh value
%    %thresh = mean([1, min(min(matin(round(.2*end):round(.8*end),
%    %round(.2*end):round(.8*end))))]); % half width
%     thresh = .1; % half width
%     
%     peaks = sqrt(col(matin < thresh).^2 + row(matin < thresh).^2);
%     
%     d = min(peaks(:));
%     
% end
% 
% if 2*d>=.8 * (s+1) % if the detected center peak diameter is more than 80% of the autocorrelogram side dim
%     
%     no_center_peak = 1;
%     
%     d = 0;
%     
% end
% 
% matin(sqrt(col.^2 + row.^2) <= d) = 0; % set center peak to zero
% 
% matout = matin;
% 
% mask = zeros(size(matin)); % build mask for plotting
% 
% mask(sqrt(col.^2 + row.^2) <= d) = 1;
% 
% end
% 
function matout = MakeSquare(matin)
% no longer pads, but instead cuts

a = diff(size(matin));

    if a<0 % if there are more rows than cols

        matout = matin(1+floor(-a/2):end-ceil(-a/2), :);
        
    elseif a>0 % if there are more cols than rows

        matout = matin(:, 1+floor(a/2):end-ceil(a/2));

    else
        matout = matin;
        
    end
    
end

function [matout, d, mask, no_center_peak, thresh] = RemoveCenterPeak(matin, r_in, thresh)
% [centercut, 2stdmajoraxis] = RemoveCenterPeak(matin);
%
% Removes center peak of autocorrelogram by finding point furthest from
% center that is half the peak (which is == 1 if using moserac).
%
% The above is the radius of center cutout. The since we are using an
% autocorr, the peak will be either round or an oblong shape. The ratio of
% the half width at the minor axis to the major axis could be the skewedness of
% fields to be added in a future version.
%
% andrew nov 10 2010

matout = [];
d = [];
mask = [];

import CMBHOME.Utils.*

r_in = 0;
thresh = .1;

s = size(matin, 1)-1; % square matrix

no_center_peak = 0;

[col, row] = meshgrid(-s/2:s/2, -s/2:s/2);

if r_in % if inner radius is given
    
    d = r_in;
    
    matin(sqrt(col.^2 + row.^2) <= d) = 0; % set center peak to zero

    matout = matin;

    mask = zeros(size(matin)); % build mask for plotting

    mask(sqrt(col.^2 + row.^2) <= d) = 1;

    return

end

if all(isnan(matin(:))), no_center_peak = 1; return; end

[~, inds] = extrema2(matin);

peak_centers = [col(inds), row(inds)];

d2 = sqrt(sum(peak_centers.^2,2));

d2 = sort(d2);

if length(d2)>1
    
    d2(1) = [];
    
    d2(d2-d2(1)>1.5*d2(1)) = [];
    
end

if ~isempty(d2)
    
    d = mean(d2(1:min([length(d2) 6])))/2;
    
end

% peaks = sqrt(col(matin >= thresh).^2 + row(matin>= thresh).^2);
% 
% d = sort(peaks(:));
% 
% d = d(find(diff(d)>1, 1, 'first')); % farthest radius of center peak

% edge_inds = edge(matin>=.1);
% 
% d = min(sqrt(col(edge_inds).^2 + row(edge_inds).^2));
% below used to work OK
        % img_prop = regionprops(matin>thresh, {'centroid', 'minoraxislength', 'majoraxislength', 'FilledArea'});
        % 
        % [~, ind] = min(sum((vertcat(img_prop(:).Centroid)-repmat((s+1)/2, length(img_prop), 2)).^2, 2)); % ind to center region
        % 
        % d = img_prop(ind).MinorAxisLength/2;

%d =  img_prop(ind).FilledArea * (img_prop(ind).MinorAxisLength/2) / ((img_prop(ind).MinorAxisLength/2)^2*pi);

if isempty(d) || 2*d>=.8 * (s+1) % if center peak is not a discrete thing, redefine thresh and look for closest non-thresh value
   %thresh = mean([1, min(min(matin(round(.2*end):round(.8*end),
   %round(.2*end):round(.8*end))))]); % half width
    thresh = .1; % half width
    
    peaks = sqrt(col(matin < thresh).^2 + row(matin < thresh).^2);
    
    d = min(peaks(:));
    
end

if 2*d>=.8 * (s+1) % if the detected center peak diameter is more than 80% of the autocorrelogram side dim
    
    no_center_peak = 1;
    
    d = 0;
    
end

matin(sqrt(col.^2 + row.^2) <= d) = 0; % set center peak to zero

matout = matin;

mask = zeros(size(matin)); % build mask for plotting

mask(sqrt(col.^2 + row.^2) <= d) = 1;

end

function [gridness, phi_v_cor, auto_corr] = MultipleEpochs(self, xdim, ydim, binside, std_smooth_kernel, rotate_inc, cel)

import CMBHOME.Utils.*

    if isempty(xdim) || isempty(ydim), % if binning dimensions are not specified, create them
    
        [x, y] = ContinuizeEpochs(self.x, self.y);

        xdim = min(x):self.spatial_scale^-1*binside:max(x); %edges of x and y dimensions

        ydim = min(y):self.spatial_scale^-1*binside:max(y);

    end
    
    gridness = nan(size(self.epoch,1), 1);
    phi_v_cor = nan(floor(180/rotate_inc)+1, size(self.epoch,1)*2);
    
    ac_dim = max([length(xdim), length(ydim)])-1;
    
    auto_corr = cell(size(self.epoch,1),1);
    
    epochs = self.epoch;
    
    for i = 1:size(epochs,1)
        
        self.epoch = epochs(i,:);
        
        [gridness(i), phi_v_cor(:,[(i-1)*2+1, (i-1)*2+2]), auto_corr{i}] = self.Gridness2(cel, 'xdim', xdim, 'ydim', ydim, 'std_smooth_kernel', std_smooth_kernel, 'rotate_inc', rotate_inc);
    
    end
        
end