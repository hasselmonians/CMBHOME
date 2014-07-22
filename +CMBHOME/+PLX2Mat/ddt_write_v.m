function [errCode] = ddt_write_v(filename, nch, npoints, freq, d)
import CMBHOME.PLX2Mat.*
% ddt_write_v(filename, nch, npoints, freq, d) Write data to a .ddt file
%
% [errCode] = ddt_write_v(filename, nch, npoints, freq, d)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%	nch - number of channels
%   npoints - number of data points per channel
%	freq - data frequency in Hz
%	d - [nch npoints] data array (in mV)
%
% OUTPUT:
%   errCode - error code: 1 for success, 0 for failure
[errCode] = mexPlex(25, filename, nch, npoints, freq, d)
