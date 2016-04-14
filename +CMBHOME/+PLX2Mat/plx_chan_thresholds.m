function [n,thresholds] = plx_chan_thresholds(filename)
import CMBHOME.PLX2Mat.*
% plx_chan_thresholds(filename): Read channel thresholds from a .plx file
%
% [n,thresholds] = plx_chan_thresholds(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%   thresholds - array of tresholds, expressed in raw A/D counts
%   n - number of channel

[n,thresholds] = mexPlex(9,filename);