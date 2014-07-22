function [n,names] = plx_adchan_names(filename)
import CMBHOME.PLX2Mat.*
% plx_adchan_names(filename): Read name for each a/d channel from a .plx file
%
% [n,names] = plx_adchan_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%   names - array of a/d channel name strings
%   n - number of channels

[n,names] = mexPlex(15,filename);