function [n,filters] = plx_chan_filters(filename)
import CMBHOME.PLX2Mat.*
% plx_chan_filters(filename): Read channel filter settings for each spike channel from a .plx file
%
% [n,filters] = plx_chan_filters(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%   filter - array of filter values (0 or 1)
%   n - number of channels

[n,filters] = mexPlex(10,filename);