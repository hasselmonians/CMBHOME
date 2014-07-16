function [n,gains] = plx_chan_gains(filename)
import CMBHOME.PLX2Mat.*
% plx_chan_gains(filename): Read channel gains from .plx file
%
% [gains] = plx_chan_gains(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%  gains - array of total gains
%   n - number of channels

[n,gains] = mexPlex(8,filename);