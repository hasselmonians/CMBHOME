function [nch, npoints, freq, d] = ddt_v(filename)
import CMBHOME.PLX2Mat.*
% ddt_v(filename) Read data from a .ddt file returning samples in mV
%
% [nch, npoints, freq, d] = ddt_v(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   nch - number of channels
%   npoints - number of data points for each channel
%   freq - A/D frequency
%   d - [nch npoints] data array (in mV)
[nch, npoints, freq, d] = mexPlex(20,filename);