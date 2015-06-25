function [cbar, clims] = smart_colorbar(clims, cmap_style)
% SMART_COLORBAR - Generates a colorbar for mapping holes in data
%
%   Finds the proper value for null data points to be plotted as white in an
%   IMAGESC plot, or similar.
%
%   [CBAR,CLIMS] = SMART_COLORBAR(CLIMS, CMAP_STYLE)
%
%   ARGUMENTS
%   * CLIMS: The desired colorbar limits
%   * CMAP_STYLE(Optional): Default jet(255). If a NX3 matrix compatible
%   with colormap, uses the colormap specified. If a string, analyzes the
%   string to find the colormap.
%
%   RETURNS
%   * CBAR: The new colobar with white background.
%   * CLIMS: The new color limits for the background. To make holes white,
%   set values to CLIMS(1)
%
% andrew 3 nov 2010
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com)   
if ~exist('cmap_style','var')
    cbar = jet(255);
else
    if ischar(cmap_style)
        eval(['cbar = ' cmap_style ';']); % initialize cbar
    else
        cbar = cmap_style;
    end
end

cbar = cat(1, [1 1 1], cbar);

nbins = size(cbar, 1)-1;

range = diff(clims)/(nbins-1)+diff(clims);

clims = [clims(2)-range, clims(2)];