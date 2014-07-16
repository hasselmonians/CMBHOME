function PipeRaster(cellin, varargin)

% Scatters the time course of a single element on ydim with the pipe '|'
% character. Existing plot must be 'held'

% PipeRaster(cellarr_in, ['color'], 3xN matrix of color values, ['ys'], 
% N array of y indeces for each row, ['ydim'], num)

% ARGUMENTS:
% cellin is the only required: Mx1 cell array of column vectors with your data.  It 
% will order them 1:M along the y axis

% color should be either an Mx3 matrix of row-wise 3 element color vectors,
% or a single color.

% accepts a cell array, and for each element, plots a new line of data

% andrew 8 oct 09
% update 2 feb 10 to varargin
% update 29 apr 10 to quit looping

p = inputParser;

p.addRequired('cellin', @iscell);

p.addParamValue('ys', 1:length(cellin), @(x) length(x)==length(cellin));
p.addParamValue('color', [0 0 0], @(x) size(x,2)==3);
p.addParamValue('Marker', '|', @ischar);
p.addParamValue('ydim',[-.2; .2], @(c) numel(c)==2); % upper and lower bounds on each integer
p.addParamValue('linewidth', 1, @isscalar);

p.parse(cellin, varargin{:});

ys = p.Results.ys;
color = p.Results.color;
Marker = p.Results.Marker;
ydim = p.Results.ydim;
linewidth = p.Results.linewidth;

ydim = ydim(:);

if size(cellin,2)>1 % used to include "size(cellin,1)>1 && "
    
    color = MakeColor(size(cellin)); % makes array of colors so that all cells have same jet color

    cellin = reshape(cellin, numel(cellin), 1);

end % assume this is coming from CMBHOME and is {epochs x cells}

n_units = length(cellin);

unit_length = cellfun(@length, cellin);

ydims = reshape(repmat(1:length(cellin), 2, 1) , 1, 2*length(cellin))';

ydims = ydims + repmat(ydim,n_units, 1); 

cell_ydims = cell(length(cellin), 1);

for i = 1:n_units
    
    cell_ydims{i} = repmat(ydims(i*2-1:i*2), 1, unit_length(i));

end

if strcmp(Marker, '.')
    
    ydims = [cell_ydims{:}]; 
    ydims = ydims(2,:)+ydim(2); 
    
    scatter(vertcat(cellin{:}), ydims, '.k', 'SizeData', 4); 
    ylim([0 max(ydims)+1]);
    return; 

end % plot dots if thats what was requested

line( repmat(vertcat(cellin{:})', 2, 1), [cell_ydims{:}], 'Color', 'k', 'linewidth', linewidth);

if numel(color)==3 % if one color make sure its black, or change it
    if ~isequal(color(:)',[0 0 0])
        h = findobj(gca, 'type', 'line');
        set(h(end-sum(unit_length)+1:end), 'Color', color);
    end    
elseif size(color, 1)==n_units
    
    h = findobj(gca, 'type', 'line');
    inds = cat(1, 0, cumsum(unit_length));
    
    for i = 1:n_units
        
        set(h(end-inds(i+1)+1:end-inds(i)), 'Color', color(i,:));

    end

end

ylim([min(ydims) max(ydims)])

end

function color = MakeColor(matsize)

color = zeros(prod(matsize), 3);

colors = jet(matsize(2))*.85;

for i = 1:matsize(2)
    
    color((i-1)*matsize(1)+1:i*matsize(1), :) = repmat(colors(i,:), matsize(1), 1);
    
end

end