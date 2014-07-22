function [nch, npoints, freq, d] = ddt(filename)
import CMBHOME.PLX2Mat.*
% ddt(filename) Read data from a .ddt file
%
% [nch, npoints, freq, d] = ddt(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%
% OUTPUT:
%   nch - number of channels
%   npoints - number of data points for each channel
%   freq - A/D frequency
%   d - [nch npoints] data array 
[nch, npoints, freq, d] = mexPlex(1,filename);