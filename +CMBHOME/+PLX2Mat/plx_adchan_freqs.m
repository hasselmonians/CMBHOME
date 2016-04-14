function [n,freqs] = plx_adchan_freqs(filename)
import CMBHOME.PLX2Mat.*
% plx_adchan_freq(filename): Read the per-channel frequencies for analog channels from a .plx file
%
% [n,freqs] = plx_adchan_freq(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%   freqs - array of frequencies
%   n - number of channels

[n,freqs] = mexPlex(12,filename);