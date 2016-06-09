function [gridness, props]  = Gridness(self, cel, varargin)
% Calculates the original gridness score of 'cel'. Continuizes all epochs.
%
% ARGUMENTS
%
%   cel             1x2 vector indicating tetrode index and cell index
%
% RETURNS
%
%   gridness        a number indicating gridness score
%   props           struct of grid cell properties (see below)
%
%   GRID CELL PROPERTY STRUCT FIELDS:
%
    %   periodicity(n, 2)       vector of the autocorrelogram correlation as a function of
    %                           rotation (see ALGORITHM below)
    %   
    %   rate_map         ratemap used in analysis
    %
    %
    %   auto_corr       matrix of the rate map autocorrelogram with donut cut
    %                   out of it. There is a second matrix along the third
    %                   dimension indicating where cuts were made. This can be
    %                   used in the AlphaData property of an image plot
    %
    %   auto_corr_mask  binary matrix where 1 indicates values used in analysis
    %                   (donut)
    %
    %   eccentricity    eccentricity of ellipse created by first six peaks in
    %                   auto_corr around center
    %
    %   e_angle         angle of rotation of major axis of ellipse (useful when
    %                   setting correction)
    %
    %   e_skew          ratio of major:minor axis of ellipse. This is the
    %                   scaling along the major axis when performing correction 
    %                   with grid3.
    %
    %   r_out           radius of outer circle (pixels)
    %
    %   r_in            radius of inner circle (pixels)
    %
    %   spacing         field spacing in pixels
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
%   grid3               (0) if 1, solves for elipse that is grid field, and
%                       then distorts it to create circle
%   r_out               radius of the outer donut
%   r_in                radius of the inner donut
%   continuize_epochs   (0) if 1, returns the grid score for rate map of
%                       cumulative data for all root.epoch(s). if 2,
%                       returns grid score for each root.epoch separately
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
% v1 dec 2009 - based upon orig gridness score (conj grid cells by sargolini et al)
% v2 may 6 2010 - added eliptical correction (grid 3 param) 
%   syntax: [gridness, phi_v_cor, auto_corr, eccentricity, grid3rotate, grid3factor, r_out, r_in]  = Gridness(self, cel, varargin)
% v3 jan 10 2011 - created grid property struct to simplify output
%
% [gridness, props] = root.Gridness(cel)
% [gridness, props] = root.Gridness(cel, params)
p = inputParser;

p.addRequired('self')
p.addRequired('cel', @isnumeric)
p.addParamValue('xdim', [], @isnumeric);
p.addParamValue('ydim', [], @isnumeric);
p.addParamValue('clims', [], @(c) numel(c)==1);
p.addParamValue('continuize_epochs', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('supress_plot', 1, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('figure_handle', [], @(c) numel(c)==1);
p.addParamValue('std_smooth_kernel', 4, @isnumeric);
p.addParamValue('binside', 3, @isnumeric)
p.addParamValue('rotate_inc', 3, @isnumeric)
p.addParamValue('grid3', 0, @isnumeric)
p.addParamValue('r_out', 0, @isnumeric)
p.addParamValue('r_in', 0, @isnumeric)
p.addParamValue('thresh', .1, @isnumeric)
p.addParamValue('grid3rotate', 0, @isnumeric)
p.addParamValue('grid3factor', 0, @isnumeric)
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
grid3 = p.Results.grid3;
r_out = p.Results.r_out;
r_in = p.Results.r_in;
thresh = p.Results.thresh;
grid3rotate = p.Results.grid3rotate;
grid3factor = p.Results.grid3factor;
autocorr = p.Results.autocorr;

import CMBHOME.Utils.* % we need the super cool extrema2 function

gridness = [];

props.periodicity = NaN; % initialize grid cell properties
props.rate_map = NaN;
props.rate_map_mask = NaN;
props.orig_auto_corr = NaN;
props.auto_corr = NaN;
props.auto_corr_mask = NaN;
props.eccentricity = NaN;
props.e_angle = NaN;
props.e_skew = NaN;
props.r_out = NaN;
props.r_in = NaN;
props.angle = NaN;
props.spacing = NaN;
props.ellipse = NaN;

if ~isempty(self)
    if size(self.epoch, 1)>1 && continuize_epochs==0 % use some recursion to do multiple epochs

        [gridness, props] = MultipleEpochs(self, xdim, ydim, binside, std_smooth_kernel,...
                                            rotate_inc, cel, grid3, r_out, r_in, thresh, grid3factor, grid3rotate);

        if ~supress_plot, PlotIt(props.auto_corr, props.periodicity, gridness); end

        return

    end
end

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

props.orig_auto_corr = auto_corr;

[auto_corr, peak_radius, mask, no_center_peak] = RemoveCenterPeak(auto_corr, r_in, thresh);

if no_center_peak, disp('No center peak'), return, end

[succ, props] = RemoveOuterCircle(auto_corr, peak_radius, mask, thresh, grid3, binside, r_out, grid3rotate, grid3factor, props);

[~, ~, props.periodicity] = AutoCorrRotation(rot90(props.auto_corr, 3), props.auto_corr, 'cut_circle', 0, 'supress_plot', 1, 'rotate_inc', rotate_inc);

gridness = min(props.periodicity(rot_ind_maxs, 2))-max(props.periodicity(rot_ind_mins, 2)); % solve for gridness

props.periodicity(:,1) = props.periodicity(:,1) + 90; % shift correlation angles to 0-180 degree

props.rate_map = rate_map;
props.rate_map_mask = rmmask;
props.r_in = peak_radius;

if ~supress_plot, PlotIt(props.auto_corr, props.periodicity, gridness); end

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

function [succ, props] = RemoveOuterCircle(auto_corr, peak_width, mask, thresh, grid3, binside, r_out, grid3rotate, grid3factor, props)
%
% auto_corr         matrix of autocorrelation
% d                 half the width of peaks (radius of peak edge)
% mask              marks where cuts in auto_corr were made
% thresh            threshold for peaks

import CMBHOME.Utils.* % for extrema2 function

eccentricity = NaN;

spacing = NaN;

succ = 0;

s = size(auto_corr, 1)-1; % square matrix

[col, row] = meshgrid(-s/2:s/2, -s/2:s/2); % indices on grid

[~, inds] = extrema2(auto_corr);

peak_centers = [col(inds), row(inds)];

d2 = sqrt(sum(peak_centers.^2,2));

peak_centers(d2<1.2*peak_width, :) = [];

d2(d2<1.2*peak_width) = [];

if r_out % if outter radius is given, find peaks within radius
    
    peak_centers(d2>r_out,:) = [];
    
    d2(d2>r_out) = [];
    
    l = length(d2);
    
    d2 = d2(1:min([l, 6]));
    
    peak_centers = peak_centers(1:min([l, 6]), :);
    
    if isempty(d2), d2 = s/2; end
        
elseif ~isempty(d2)
    
    succ = 1;

    [d2, ind] = sort(d2(:));
    
    peak_centers = peak_centers(ind, :);
    
    peak_centers(d2-d2(1)>.6*d2(1),:) = []; % clear peaks much farther away from the first
    
    d2(d2-d2(1)>.6*d2(1)) = [];
     
    l = length(d2);
    
    d2 = d2(1:min([l, 6]));
    
    peak_centers = peak_centers(1:min([l, 6]), :);
    
    if isempty(d2), d2 = s/2; end
    
else
    
    disp('No surrounding peaks.');
    
    d2 = s/2;
    
end

if (grid3 & length(d2)==6) || (grid3 & r_out & grid3rotate & grid3factor) 
    
    [auto_corr, grid3factor, mask, eccentricity, grid3rotate, cs, r_out] = Ovalate(auto_corr, d2, peak_width, peak_centers, mask, r_out, grid3rotate, grid3factor);
    
    props.ellipse = cs; % assign ellipse coefficients
    props.r_out = r_out;
    
elseif r_out % user spec outer radius but 6 peaks were not detected, or
      
    [auto_corr, mask ] = MakeDonutAndMask(auto_corr, mask, r_out);
    props.r_out = r_out;
        
else

    [auto_corr, mask ] = MakeDonutAndMask(auto_corr, mask, d2(end)+peak_width);
    props.r_out = d2(end)+peak_width;
    
end

props.auto_corr = auto_corr;
props.auto_corr_mask = ~mask;
props.eccentricity = eccentricity;
props.e_angle = grid3rotate;
props.e_skew = grid3factor;

end

function [auto_corr, grid3factor, mask, eccentricity, grid3rotate, cs, r_out] = Ovalate(auto_corr, d2, peak_width, peak_centers, mask, r_out, grid3rotate, grid3factor)
% This function takes the 6 peaks in the autocorrelogram, fits an adapted
% oval equation to it, and then squishes the autocorrelation to make the
% pattern roouunund

%figure('position', [0 0 1300 600])

import CMBHOME.Utils.*

rot_ind_mins = round([30/30, 90/30, 150/30])+1; % indexes in rotated ac to check

rot_ind_maxs = round([60/30, 120/30])+1;

s = size(auto_corr, 1);

if r_out & grid3rotate & grid3factor % if we are to assume an elliptical correction
    
    auto_corr = imrotate(auto_corr, grid3rotate, 'bicubic', 'crop'); % rotate auto_corr so that major axis is up
    mask = imrotate(mask, grid3rotate - 90, 'bicubic', 'crop');
    auto_corr = imresize(auto_corr, [ceil(s*grid3factor) s], 'bicubic'); % resize it
    mask = imresize(mask, [ceil(s*grid3factor) s], 'bicubic');
    
    [auto_corr, mask] = MakeDonutAndMask(auto_corr, mask, r_out, peak_width);
    
    eccentricity = sqrt(1 - factor^2);
    
else
    
    % Try without correction
    
    [nocorrect.auto_corr, nocorrect.mask] = MakeDonutAndMask(auto_corr, mask, min([d2(end)+peak_width/2 s/2]));
    [~, ~, nocorrect.phi_v_cor] = AutoCorrRotation(rot90(nocorrect.auto_corr, 3), nocorrect.auto_corr, 'cut_circle', 0, 'supress_plot', 1, 'rotate_inc', 30);
    nocorrect.gridness = min(nocorrect.phi_v_cor(rot_ind_maxs, 2))-max(nocorrect.phi_v_cor(rot_ind_mins, 2));
    
    % Fit ellipse and solve for minor and major axis, etc

    cs = EllipseDirectFit(peak_centers); % coefficients for equation for an ellipse

    [A, B, phi, x, y] = EllipseFromCoef(cs, peak_centers);
    
    ec.factor = min(A, B)/max(A,B);
    
    if ec.factor >.5 % if major axis is less than twice the length of the minor axis
    
        ec.auto_corr = imrotate(auto_corr, rad2deg(phi)+90, 'bicubic', 'crop'); % rotate auto_corr so that major axis is up
        ec.mask = imrotate(mask, 90+rad2deg(phi), 'bicubic', 'crop');
        ec.auto_corr = imresize(ec.auto_corr, [ceil(s*ec.factor) s], 'bicubic'); % resize it
        ec.mask = imresize(ec.mask, [ceil(s*ec.factor) s], 'bicubic');
        
        [ec.auto_corr, ec.mask ] = MakeDonutAndMask(ec.auto_corr, ec.mask, B+peak_width/2, peak_width);
        [~, ~, ec.phi_v_cor] = AutoCorrRotation(rot90(ec.auto_corr, 3), ec.auto_corr, 'cut_circle', 0, 'supress_plot', 1, 'rotate_inc', 30);
        ec.gridness = min(ec.phi_v_cor(rot_ind_maxs, 2))-max(ec.phi_v_cor(rot_ind_mins, 2));
        ec.eccentricity = sqrt(1-(B/A)^2);
        
    else
        
        ec.gridness = -100; % some arbitrarily low number
        ec.auto_corr = auto_corr;
        ec.mask = mask;
        cs = NaN;
        
    end
   
%     tmph = gcf;
%         
%     figure, subplot(1, 2, 1), 
%     imagesc(auto_corr), axis equal, set(gca, 'ydir', 'normal'), hold on
%     scatter(peak_centers(:,1)+s/2+.5, peak_centers(:,2)+s/2+.5, 'g*'),
%     line([s/2+.5, s/2+A*cos(phi)+.5], [s/2+.5, s/2+A*sin(phi)+.5], 'Color', 'r', 'LineWidth', 2);
%     line([s/2+.5, s/2+B*cos(phi+pi/2)+.5], [s/2+.5, s/2+B*sin(phi+pi/2)+.5], 'Color', 'b', 'LineWidth', 2); title('major app'), axis equal, hold off
%     line(x+s/2+.5, y+s/2+.5, 'color', 'k', 'linewidth', 1.5) 
%     
%     subplot(1, 2, 2), imagesc(ec.auto_corr), title(['g: ' num2str(ec.gridness)]), set(gca, 'ydir', 'normal'), axis equal
% 
%     %print(gcf, '-dpsc2', 'grid3.ps', '-append')
%     %close(gcf)
%     figure(tmph)
     
    if ec.gridness < nocorrect.gridness
    
        grid3rotate = 90; % no rotation
        grid3factor = 1;
        auto_corr = nocorrect.auto_corr;
        mask = nocorrect.mask;
        eccentricity = 0;
        r_out = min([d2(end)+peak_width/2 s/2]);
        
    else
        
        grid3rotate = phi;
        grid3factor = ec.factor;
        auto_corr = ec.auto_corr;
        mask = ec.mask;
        eccentricity = ec.eccentricity;
        r_out = B+peak_width/2;
        
    end
    
end

end

function [matout, d, mask, no_center_peak] = RemoveCenterPeak(matin, r_in, thresh)
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

function [auto_corr, mask] = MakeDonutAndMask(auto_corr, mask, r, recut_center)

    if ~exist('recut_center', 'var'), recut_center = 0; end
    
    s1 = size(auto_corr, 1)-1;
    s2 = size(auto_corr, 2)-1;

    [col, row] = meshgrid(-s2/2:s2/2, -s1/2:s1/2); % indices on grid % cut the donut

    r = min([ (s1-1)/2, r]); % set outer radius

    auto_corr(sqrt(col.^2 + row.^2) > r) = 0;
    
    if recut_center
        auto_corr(sqrt(col.^2 + row.^2)<=recut_center) = 0;

        mask(sqrt(col.^2 + row.^2)<=recut_center) = 1;
    end
    
    mask(sqrt(col.^2 + row.^2) > r) = 1;
    
    r = ceil(r); % reshape to square
    
    ds = size(auto_corr, 1)/2-r;
    
    ds2 = size(auto_corr, 2)/2-r;
    
    if ds<1, ds = .9; end
    if ds2<1, ds2 = .9; end
    
    auto_corr = auto_corr(ceil(ds):end-floor(ds),ceil(ds2):end-floor(ds2));
    
    mask = mask(ceil(ds):end-floor(ds),ceil(ds2):end-floor(ds2));
    
    auto_corr = MakeSquare(auto_corr);
    
    mask = MakeSquare(mask);
    
    mask = round(mask);
    
end

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

function [gridness, props] = MultipleEpochs(self, xdim, ydim,...
                                            binside, std_smooth_kernel, rotate_inc,...
                                            cel, grid3, r_out, r_in, thresh,...
                                            grid3factor, grid3rotate)

    import CMBHOME.Utils.*

    if isempty(xdim) || isempty(ydim), % if binning dimensions are not specified, create them
    
        [x, y] = ContinuizeEpochs(self.x, self.y);

        xdim = min(x):self.spatial_scale^-1*binside:max(x); %edges of x and y dimensions

        ydim = min(y):self.spatial_scale^-1*binside:max(y);

    end
    
    gridness = nan(size(self.epoch,1), 1);
    
    epochs = self.epoch;
    
    for i = 1:size(epochs,1)
        
        self.epoch = epochs(i,:);

        [gridness(i), props(i)] = self.Gridness(cel, 'xdim', xdim, 'ydim', ydim, 'std_smooth_kernel', std_smooth_kernel, 'rotate_inc', rotate_inc, 'grid3', grid3', 'r_out', r_out, 'r_in', r_in, 'thresh', thresh, 'grid3factor', grid3factor, 'grid3rotate', grid3rotate);
    
    end
        
end