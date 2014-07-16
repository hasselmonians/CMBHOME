function plot_trajectory(self, cel, varargin)
% root.plot_trajectory(cel)
%
% Plots animal trajectory for root.epochs (if multiple, they are all
% plotted on top of one another), along with spikes overlayed in red dots.
%
% andrew 3 april 2010

p = inputParser;

p.addRequired('self')
p.addRequired('cel', @isnumeric)
p.addParamValue('show_scale', 1, @isnumeric)
p.addParamValue('ScatterParams', {'.r', 'SizeData', 25}, @iscell);

p.parse(self, cel, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
show_scale = p.Results.show_scale;
ScatterParams = p.Results.ScatterParams;
    
    import CMBHOME.Utils.*

    pad = [-.03 .02]; % percent padding around plot
    
    [spk_x, spk_y] = ContinuizeEpochs(self.spk_x(cel), self.spk_y(cel));
        
    x = Cell2MatPlus(self.x);
    y = Cell2MatPlus(self.y);
    
    line(x, y, 'Color', 'k'), hold on;
       
    scatter(spk_x, spk_y, ScatterParams{:}), hold off

    axis equal
    
    axis off
    
    set(gca, 'Box', 'on')

    xs = [min(min(x)) max(max(x))];
    ys = [min(min(y)) max(max(y))];
    
    xlim(diff(xs).*pad+xs);
    ylim(diff(ys).*pad+ys);
    
    if show_scale
        line([xs(1)+.75*diff(xs), xs(2)], [-.03*diff(ys)+ys(1), -.03*diff(ys)+ys(1)], 'Color', 'k', 'LineWidth', 2);
        text(xs(2), -.03*diff(ys)+ys(1), [num2str(.25*diff(xs)*self.spatial_scale, 3) ' cm'], 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom');
    end
    
    %title('Animal Trajectory with Spikes');

end
