function [n, npw, ts, wave] = plx_waves(filename, ch, u)
import CMBHOME.PLX2Mat.*
% plx_waves(filename, channel, unit): Read waveform data from a .plx file
%
% [n, npw, ts, wave] = plx_waves(filename, channel, unit)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based channel number
%   unit  - unit number (0- unsorted, 1-4 units a-d)
% OUTPUT:
%   n - number of waveforms
%   npw - number of points in each waveform
%   ts - array of timestamps (in seconds) 
%   wave - array of waveforms [npw, n], raw a/d values

[n, npw, ts, wave] = mexPlex(6,filename, ch, u);
