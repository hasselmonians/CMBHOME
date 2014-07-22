function [n,gains] = plx_adchan_gains(filename)
import CMBHOME.PLX2Mat.*
% plx_adchan_gains(filename): Read analog channel gains from .plx file
%
% [n,gains] = plx_adchan_gains(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%  gains - array of total gains
%  n - number of channels

[n,gains] = mexPlex(11,filename);