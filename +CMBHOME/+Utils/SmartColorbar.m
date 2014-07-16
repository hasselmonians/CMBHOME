function [cbar, clims] = SmartColorbar(clims, cmap_style)
% CMBHOME.Utils.SmartColorbar(clims, cmap_style)
%
% finds the proper value for null data points to be plotted as white in an
% IMAGESC plot, or similar.
%
% ARGUMENTS
%   clims           (1x2 vector) either the range of the data, or the prefered range of the
%                   colorbar
%   cmap_style      string. ex. 'jet(255)'
%
% RETURNS
%   cbar            Nx3 array of colorbar values. lowest is white
%   clims           new clims. clims(1) should be included as the value of
%                   the null data points to be properly mapped to white
%
% andrew 3 nov 2010

eval(['cbar = ' cmap_style ';']); % initialize cbar

cbar = cat(1, [1 1 1], cbar);

nbins = size(cbar, 1)-1;

range = diff(clims)/(nbins-1)+diff(clims);

clims = [clims(2)-range, clims(2)];