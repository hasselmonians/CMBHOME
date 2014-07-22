function [n,samplecounts] = plx_adchan_samplecounts(filename)
import CMBHOME.PLX2Mat.*
% plx_adchan_samplecounts(filename): Read the per-channel sample counts for analog channels from a .plx file
%
% [n,samplecounts] = plx_adchan_samplecounts(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%   n - number of channels
%   samplecounts - array of sample counts

[n,samplecounts] = mexPlex(23,filename);