function [phi, max_cor, phi_v_cor] = AutoCorrRotation(ac1, ac2, varargin)
% [phi] = AutoCorrRotation(ac1, ac2);
% [phi, max_cor, phi_v_cor] = AutoCorrRotation(ac1, ac2);
%
% Calculates angle phi for which ac1 and ac2 are maximally correlated.
% Useful for checking for rotation of grid field orientation between
% conditions
%
% ARGUMENTS
%   ac1         matrix of autocorrelogram 1
%   ac2         matric of autocorrelogram 2
%
% RETURNS
%   phi         angle of max correlation (degrees)
%   max_cor     value of max correlation (0-1)
%   phi_v_cor   Nx2 vector [phis; cors]
%
% PARAMETERS
%   cut_circle          0 or 1 (1), cuts a circle out from each ac so that corners
%                       do not affect rotated score
%   rotate_inc          (default 1 degree). what increment for rotation?
%   supress_plot        0 or 1 (0). plots a subplot with both
%                       autocorrelations to be compared (after cut, if performed), as well as
%                       periodicity of correlation
%
% ALGORITHM
%   1. Unless cut_circle==0, finds smallest dimenstion between both ac's,
%   and cuts a circle out from each ac, to a matching diameter. If
%   cut_circle==0, ac1 and ac2 must be the same size
%
%   2. calculates the pearsons correlation coefficient for rotations of 1
%   degree counterclockwise of ac2 against ac1 (from -90 to 90, which makes
%   sense only if both ac1 and ac2 are autocorrelations
%
%   3. Finds max of the correlation score, and returns both vectors, as
%   well.
%   

import CMBHOME.Utils.*

p = inputParser;

p.addRequired('ac1', @isnumeric)
p.addRequired('ac2', @isnumeric)
p.addParamValue('cut_circle', 1, @(c) c==1 | c==0);
p.addParamValue('rotate_inc', 1, @isnumeric)
p.addParamValue('supress_plot', 0, @(c) c==1 | c==0);

p.parse(ac1, ac2, varargin{:});

ac1 = p.Results.ac1;
ac2 = p.Results.ac2;
cut_circle = p.Results.cut_circle;
rotate_inc = p.Results.rotate_inc;
supress_plot = p.Results.supress_plot;

if cut_circle, [ac1, ac2] = CutCircle(ac1, ac2); end

if size(ac1)~=size(ac2), error('AutoCorrRotation: ac1 and ac2 must be the same size'); end

p_ac1 = ac1; % for plots at the end
p_ac2 = ac2;

ac1(isnan(ac1)) = 0;
ac2(isnan(ac2)) = 0;

phi = -90:rotate_inc:90;

phi_v_cor = zeros(length(phi), 2);

for i = 1:length(phi)
    
    rot_ac = imrotate(ac2, phi(i), 'bilinear', 'crop');
    
    r = corrcoef(rot_ac(:), ac1);

    phi_v_cor(i, 2) = r(2);
                
end

phi_v_cor(:,1) = phi;

[max_cor, ind] = max(phi_v_cor(:,2));

phi = phi(ind);

if ~supress_plot, PlotIt(phi_v_cor, p_ac1, p_ac2, phi, max_cor); end

end

function PlotIt(phi_v_cor, ac1, ac2, phi, max_cor)

    figure('Position', [50 50 800 500]);

    subplot(2, 2, 1), imagesc(ac1, 'AlphaData', ~isnan(ac1)), axis square, axis off
    subplot(2, 2, 2), imagesc(ac2, 'AlphaData', ~isnan(ac2)), axis square, axis off

    subplot(2, 2, 3:4)
    line(phi_v_cor(:,1), phi_v_cor(:,2), 'Color', 'k', 'LineWidth', 2.5), hold on
    scatter(phi, max_cor, 'or', 'SizeData', 45), hold on
    text(phi+2, max_cor, ['Phi Max Correlation: ' num2str(phi)]); 
    xlim([-90 90]);

end

function [ac1_cut, ac2_cut] = CutCircle(ac1, ac2)

    min_d = min([size(ac1), size(ac2)]);
    
    ac1_cut = nan(min_d);
    ac2_cut = nan(min_d);
    
    ac1 = ac1(  floor((size(ac1, 1)-min_d)/2)+1:end-ceil((size(ac1, 1)-min_d)/2),...
                floor((size(ac1, 2)-min_d)/2)+1:end-ceil((size(ac1, 2)-min_d)/2));

    ac2 = ac2(  floor((size(ac2, 1)-min_d)/2)+1:end-ceil((size(ac2, 1)-min_d)/2),...
                floor((size(ac2, 2)-min_d)/2)+1:end-ceil((size(ac2, 2)-min_d)/2));

    [c, r] = meshgrid(  -(size(ac1, 2)-1)/2:(size(ac1, 2)-1)/2,... 
                        -(size(ac1, 1)-1)/2:(size(ac1, 1)-1)/2 );
    
    ac1_cut(sqrt(c.^2 + r.^2) < min_d/2) = ac1(sqrt(c.^2 + r.^2) < min_d/2); % cut ac1
    
    ac2_cut(sqrt(c.^2 + r.^2) < min_d/2) = ac2(sqrt(c.^2 + r.^2) < min_d/2); % cut ac2
    
end