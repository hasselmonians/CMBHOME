function [rate_map, xs, ys, xdim, ydim, no_occupancy] = plot_rate_map(self, cel, suppress_plot, xdim, ydim, clims, mergeepochs)
% [rate_map, xs, ys] = root.plot_rate_map(cel, suppress_plot, clims)
%
% Plots spatial rate map for cell cel = [ tetrode, cell ] for all root.epoch
% (multiple epochs concatinated).
%
% Returns matrix rate_map with x and y dimensions in vectors xs and ys. If
% suppress_plot = 1, the rate map is not plotted, but just calculated
% 
% andrew bogaard 3 april 2010


    if ~exist('suppress_plot', 'var')
        suppress_plot=0;
    end 
%   xdim                vector of bin edges along x dimension
%   ydim                vector of bin edges along y dimension
    
    if ~exist('mergeepochs', 'var'), mergeepochs = 1; end
    
    if ~exist('xdim', 'var') || ~ exist('ydim', 'var')
        xdim = [];
        ydim = [];
    end
    
    pad = [-.05 .05]; % percent pad plot

    import CMBHOME.Utils.*

    [occupancy, xdim, ydim] = self.Occupancy(xdim, ydim, mergeepochs);

    no_occupancy = occupancy~=0; % mark indeces where there was no occupancy so we can correct after smoothing
    
    if mergeepochs
        [spk_x, spk_y] = ContinuizeEpochs(self.cel_x, self.cel_y); 
    else
        spk_x = self.cel_x;
        spk_y = self.cel_y;
    end

    if ~iscell(spk_x)
        spikes = hist3([spk_x, spk_y], 'Edges', {xdim, ydim});

        rate_map = SmoothMat(spikes, [5, 5], 1)./SmoothMat(occupancy, [5, 5], 1); % smooth the spikes and occupancy with a 5x5 bin gaussian with std=1
        
        rate_map(~no_occupancy) = 0; % set no occupancy to zero
    
        rate_map = rate_map'; % returns these three
        xs = [min(xdim) max(xdim)];
        ys = [min(ydim) max(ydim)];
        no_occupancy = no_occupancy';
        
    else % multiple epochs
        
        rate_map = zeros(size(occupancy, 2), size(occupancy, 1), size(occupancy, 3));
        
        new_occupancy = zeros(size(occupancy, 2), size(occupancy,1), size(occupancy,3));
        
        for i = 1:size(occupancy,3)
        
            spikes = hist3([spk_x{i}, spk_y{i}], 'Edges', {xdim, ydim});

            tmp = SmoothMat(spikes', [5, 5], 1)./SmoothMat(occupancy(:,:,i)', [5, 5], 1); % smooth the spikes and occupancy with a 5x5 bin gaussian with std=1

            tmp(~no_occupancy(:,:,i)') = 0;
            
            rate_map(:, :, i) = tmp;

            xs = [min(xdim) max(xdim)];
            ys = [min(ydim) max(ydim)];
            new_occupancy(:,:,i) = no_occupancy(:,:,i)';
            
        end
        
    end
        
    
    if (~exist('clims', 'var') || isempty(clims))
        clims = [0 max(rate_map(:))];
    end
    
    if ~suppress_plot && mergeepochs==1
        figure
        t=imagesc(xdim, ydim, rate_map, clims);
        colormap jet(255);
        axis equal

        set(gca, 'Box', 'on')

        xlim(xs.*pad+xs);
        ylim(ys.*pad+ys);

        text(xs(2)+.07*xs(2), 1.02*ys(2), [int2str(max(max(rate_map))) 'Hz'], 'FontSize', 12, 'FontWeight', 'bold');
        
        rmpos=get(gca, 'Position');

        set(gca,'YDir','normal'); % so plotting functions dont reverse axis

        set(t, 'AlphaDataMapping', 'none');   %   until this works properly
                                                  %without screwing up other axes in the subplot, screw it!
        set(gca, 'DrawMode', 'fast');
        set(t,'AlphaData', no_occupancy);

        set(gca,'XTickLabel',cellfun(@(x) str2num(x), get(gca,'XTickLabel'))*(self.spatial_scale))
        set(gca,'YTickLabel',cellfun(@(x) str2num(x), get(gca,'YTickLabel'))*(self.spatial_scale))
        
    elseif ~suppress_plot && mergeepochs == 0
        for i = 1:size(occupancy,3)
            figure
            t=imagesc(xdim, ydim, rate_map(:,:,i), clims);
            colormap jet(255);
            axis equal

            set(gca, 'Box', 'on')

            xlim(xs.*pad+xs);
            ylim(ys.*pad+ys);

            text(xs(2)+.07*xs(2), 1.02*ys(2), [int2str(max(max(rate_map))) 'Hz'], 'FontSize', 12, 'FontWeight', 'bold');

            rmpos=get(gca, 'Position');

            set(gca,'YDir','normal'); % so plotting functions dont reverse axis

            set(t, 'AlphaDataMapping', 'none');   %   until this works properly
                                                      %without screwing up other axes in the subplot, screw it!
            set(gca, 'DrawMode', 'fast');
            set(t,'AlphaData', no_occupancy(:,:,i)');
            
            %set(gca,'XTickLabel',str2num(get(gca,'XTickLabel'))*(self.spatial_scale))
            %set(gca,'YTickLabel',str2num(get(gca,'YTickLabel'))*(self.spatial_scale))
        end
    end

end
