function [n,names] = plx_chan_names(filename)
import CMBHOME.PLX2Mat.*
% plx_chan_names(filename): Read name for each spike channel from a .plx file
%
% [n,names] = plx_chan_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%   names - array of channel name strings
%   n - number of channels

[n,names] = mexPlex(14,filename);