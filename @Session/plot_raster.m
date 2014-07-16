function plot_raster(self, cel, varargin)
%
% Peristimulus time histogram of spikes between times
%
% (1) root.plot_raster(cel);
% (2) root.plot_raster(cel, epochs_offsets);
% (3) root.plot_raster(cel, epochs_offsets, plot_pd)
%
% (1) Prints to current axes raster plot of spikes in cel ([tetrode index, cell
% index]), between epochs in root.epoch, and lines up all epochs by their start time
% (2) Prints to current axes raster plot of spikes in cel between epochs in
% root.epoch, and lines them up by epochs_offsets, which is a N element 
% vector as long as root.epoch. This allows you to center justify, event
% justify, etc. 
% (3) If plot_pd has a non-zero value, a probability density function will be plotted
% below the raster at binsize plot_histogram
%
% In laymens: the function without epochs_offsets assumes that the 0 time
% for each epoch is the starting time of the epoch. If you include
% epochs_offsets, the 0 time for each epoch is defined in epochs_offsets

import CMBHOME.Utils.*

p = inputParser;

p.addRequired('self');
p.addRequired('cel');
p.addOptional('epochs_offsets', [], @(x) isnumeric(x));
p.addOptional('plot_pd',  0, @(x) numel(x)==1);
p.addParamValue('Justify', 'left', @(x) any(strcmp(x, {'left', 'right', 'center'})));
p.addParamValue('Marker', '|', @ischar);
p.addParamValue('Axis', {''}, @iscellstr);
p.addParamValue('CenterLine', 0, @(x) numel(x)==1);
p.addParamValue('Color', [0 0 0], @(x) size(x,2)==3); 
p.addParamValue('ydim',[-.2; .2], @(c) numel(c)==2); % upper and lower bounds on each integer

p.parse(self, cel, varargin{:});

epochs_offsets = p.Results.epochs_offsets;
plot_pd = p.Results.plot_pd;
Justify = p.Results.Justify;
Marker = p.Results.Marker;
varAxis = p.Results.Axis;
CenterLine = p.Results.CenterLine;
Color = p.Results.Color;
ydim = p.Results.ydim;

if isempty(epochs_offsets)
    switch Justify
        case 'left'
            epochs_offsets = self.epoch(:,1); % if no offsets, offset by starting times
        case 'center'
            epochs_offsets = (self.epoch(:,2)-self.epoch(:,1))/2 + self.epoch(:,1);
        case 'right'
            epochs_offsets = self.epoch(:,2); % offset by ending times
    end
end

% build cell array of firing

cell_spk = self.spk_ts(cel); % cell array of spike times

if iscell(cell_spk)
    for i = 1:size(cell_spk,1) % take care of multiple epochs
        for j = 1:size(cell_spk,2) % take care of multiple cells
            cell_spk{i, j} = cell_spk{i, j} - epochs_offsets(i);
        end
    end
else

    cell_spk = {cell_spk};

end

if plot_pd
    
    pos = get(gca, 'Position');
    
    set(gca,'Position', [pos(1) pos(2)+.25*pos(4) pos(3) pos(4)*.8]);
    
    PipeRaster(cell_spk, 'Marker', Marker, 'Color', Color, 'ydim', ydim);
    
    DrawCenter(CenterLine); % if indicated, draw center line
    
    SetAxis(varAxis);
        
    cell_spk = vertcat(cell_spk{:}); % make into vector
    xs = xlim;
    count = histc(cell_spk, min(cell_spk):plot_pd:max(cell_spk));
    count(end)=[];
    count = count./sum(count);
    
    axes('Position', [pos(1), pos(2)-.05*pos(4), pos(3), pos(4)*.2]);

    plot(count, 'Color', 'r', 'LineWidth', 1.5);

    if length(count)>1, xlim([1, length(count)]); end
    ylim([0 max([.15 max(count)])])
    axis off
    
else
    
    PipeRaster(cell_spk, 'Marker', Marker, 'Color', Color, 'ydim', ydim);
    
    DrawCenter(CenterLine); % if indicated, draw center line
    
    SetAxis(varAxis);
    
end
        
end

function SetAxis(varAxis)

for i=1:length(varAxis)
switch varAxis{i}
    
    case 'off'
        
        axis off
        
    case 'square'
        
        axis square
        
end
end
end

function DrawCenter(CenterLine)

if CenterLine
    
    ys = ylim;
    xs = xlim;
    
    hold on % hold axis on
    
    line([mean(xs) mean(xs)], ys, 'Color', 'r', 'LineWidth', 1, 'LineStyle', ':');
    
end
end
