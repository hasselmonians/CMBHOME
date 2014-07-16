function [rate_map, xdim, ydim, occupancy, no_occupancy] = RateMap(self, varargin)
% [rate_map, xs, ys] = root.plot_rate_map(cel, varargin)
%
% Plots spatial rate map for cell cel = [ tetrode, cell ]
%
% Returns matrix rate_map with x and y dimensions in vectors xs and ys. 
%
% ARGUMENTS
%   cel             1 x 2 vector like [tetrode index, cell index]
%   params          (see below)
%
% RETURNS
%   rate_map        matrix of F/bin (if continuize_epochs == 0), a third
%                   dimension exists for each epoch
%   xdim            vector of x dimensions
%   ydim            vector of y dimensions
%   occupancy       matrix the size of rate_map indicating occupancy in
%                   seconds
%   no_occupancy    logical matrix the size of rate_map. 1's indicate no
%                   occupancy. useful for AlphaData property of imagesc
%
% OPTIONAL PARAMETERS
%   
%   xdim                vector of bin edges along x dimension
%   ydim                vector of bin edges along y dimension
%   continuize_epochs   0 or 1 (0). If 0, ratemap is calculated for each
%                       epoch (adds 3rd dim to rate_map output, if 1, ratemap
%                       is calculated across all epochs
%   supress_plot        0 or 1 (1). If 0, plots ratemap
%   figure_handle       If supplied, plot the rate map to the given figure
%   std_smooth_kernel   STD (cm) of the gaussian kernel to smooth the rate map  
%   binside             The length in cm of the side of a bin when
%                       calculating the rate map
%   show_scale          (1) if plotting, shows scale or not
%   show_peak_f         (1) if 1, shows peak f or not
%   show_mean_f         (1) if 1, shows mean f or not
%   omit_islands        (1) if 1, solves for center of mass of occupancy,
%                       and then looks for any stray pixels and sets their 
%                       spikes to zero  
%   omit_noocc          (1) if 1, sets to zero pixels from ratemap which had
%                       occupancy values = 0
% 
% andrew bogaard 3 april 2010
% rev 2 andrew bogaard 15 sept 2010, better return variables and options

p = inputParser;

p.addRequired('self')
p.addOptional('cel', NaN, @isnumeric)
p.addParamValue('xdim', [], @isnumeric);
p.addParamValue('ydim', [], @isnumeric);
p.addParamValue('clims', [], @(c) numel(c)==2);
p.addParamValue('continuize_epochs', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('supress_plot', 1, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('figure_handle', '', @(c) numel(c)==1);
p.addParamValue('std_smooth_kernel', 3, @isnumeric);
p.addParamValue('binside', 3, @isnumeric)
p.addParamValue('show_scale', 1, @isnumeric)
p.addParamValue('show_peak_f', 1, @isnumeric)
p.addParamValue('show_mean_f', 1, @isnumeric)
p.addParamValue('omit_islands', 1, @isnumeric)
p.addParamValue('omit_noocc', 1, @isnumeric)

p.parse(self, varargin{:});

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
show_scale = p.Results.show_scale;
show_peak_f = p.Results.show_peak_f;
show_mean_f = p.Results.show_mean_f;
omit_islands = p.Results.omit_islands;
omit_noocc = p.Results.omit_noocc;
    
import CMBHOME.Utils.*

if isnan(cel), cel = self.cel; end

if any(size(cel)~=[1, 2]), error('cell must be size 1x2 like [tetrode, cell]'); end

if continuize_epochs, self = MergeEpochs(self); end % so that we don't double count data

[occupancy, xdim, ydim] = self.Occupancy(xdim, ydim, continuize_epochs, binside);

binside = self.spatial_scale * diff(xdim(1:2)); % solve for binside in case xdim and ydim were specified

if omit_islands, occupancy = OmitIslands(occupancy); end

no_occupancy = occupancy==0; % mark indeces where there was no occupancy so we can correct after smoothing
    
if continuize_epochs
    [spk_x, spk_y] = ContinuizeEpochs(self.spk_x(cel), self.spk_y(cel)); 
else
    spk_x = self.spk_x(cel);
    spk_y = self.spk_y(cel);
end

if ~iscell(spk_x)
    spikes = hist3([spk_x, spk_y], 'Edges', {xdim, ydim});

    rate_map = SmoothMat(spikes, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside)./SmoothMat(occupancy, [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside); % smooth the spikes and occupancy with a 5x5 bin gaussian with std=1

    if omit_noocc
        rate_map(no_occupancy) = 0; % set no occupancy to zero
    end

    rate_map = rate_map'; % returns these three
    
    occupancy = occupancy';

    no_occupancy = no_occupancy';
    
    mean_f = length(spk_x)/sum(self.epoch(:,2) - self.epoch(:,1));

else % multiple epochs

    rate_map = zeros(size(occupancy, 2), size(occupancy, 1), size(occupancy, 3));

    new_occupancy = zeros(size(occupancy, 2), size(occupancy,1), size(occupancy,3));

    for i = 1:size(occupancy,3)

        spikes = hist3([spk_x{i}, spk_y{i}], 'Edges', {xdim, ydim});

        tmp = SmoothMat(spikes', [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside)./SmoothMat(occupancy(:,:,i)', [5*std_smooth_kernel/binside, 5*std_smooth_kernel/binside], std_smooth_kernel/binside); 

        if omit_noocc
            tmp(no_occupancy(:,:,i)') = 0;
        end

        rate_map(:, :, i) = tmp;

        new_occupancy(:,:,i) = no_occupancy(:,:,i)';

    end

    occupancy = permute(occupancy,[2 1 3]);
    
    no_occupancy = new_occupancy;

end
        
if ~supress_plot && size(rate_map,3)==1
    PlotIt(rate_map, no_occupancy, clims, xdim, ydim, figure_handle, self, show_scale, show_peak_f, show_mean_f, mean_f);
end
    
end

function occupancy = OmitIslands(occupancy)
% Takes matrix of occupancy values, calculates center of mass of pixels>0,
% and then finds all disconnected pixels and sets them = 0

%subplot(1, 2, 1), imagesc(occupancy)

s = regionprops(occupancy>0, {'FilledArea', 'PixelIdxList'});

l = numel(s);

areas = vertcat(s.FilledArea);

[~, inds] = sort(areas);

for i = 1:length(inds)-1
    
    occupancy(s(inds(i)).PixelIdxList) = 0;
    
end 

%subplot(1, 2, 2), imagesc(occupancy)

end

function PlotIt(rate_map, no_occupancy, clims, xdim, ydim, handle, self, show_scale, show_peak_f, show_mean_f, mean_f)

    if all(isnan(rate_map)) | all(rate_map==0)
        
        text(.05, .4, 'No figure'), axis off
        
        return
        
    end

    xs = [min(xdim) max(xdim)];
    ys = [min(ydim) max(ydim)];

    if ~isempty(handle), figure(handle); end % make current the figure passed

    if isempty(clims)
        clims = [0 max(max(rate_map))];
    end
    
    pad = [-.03 .02]; % percent pad plot
    
    [cbar, clims] = CMBHOME.Utils.SmartColorbar(clims, 'jet(255)');

    rate_map(no_occupancy) = clims(1);

    imagesc(xdim, ydim, rate_map, clims); hold on;
    
    colormap(cbar);
    
    axis equal off

    xlim(diff(xs).*pad+xs);
    ylim(diff(ys).*pad+ys);
    
    if show_scale
        line([xs(1)+.75*diff(xs), xs(2)], [-.03*diff(ys)+ys(1), -.03*diff(ys)+ys(1)], 'Color', 'k', 'LineWidth', 2);
        text(xs(2), -.03*diff(ys)+ys(1), [num2str(.25*diff(xs)*self.spatial_scale, 3) ' cm'], 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom');
    end
    
    str_f = '';
    
    if show_peak_f, str_f = ['p: ' num2str(max(max(rate_map)), 2) 'Hz']; end
    
    if show_mean_f, str_f = [str_f ' m: ' num2str(mean_f, 2) 'Hz']; end
    
    if show_peak_f | show_mean_f
        text(xs(2), .045*diff(ys)+ys(2), str_f, 'FontSize', 6.8, 'FontWeight', 'bold', 'HorizontalAlign', 'right');
    end
        
    rmpos=get(gca, 'Position');

    set(gca,'YDir','normal'); % so plotting functions dont reverse axis

end