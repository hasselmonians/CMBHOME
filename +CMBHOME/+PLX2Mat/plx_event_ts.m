function [n, ts, sv] = plx_event_ts(filename, ch)
import CMBHOME.PLX2Mat.*
% plx_event_ts(filename, channel) Read event timestamps from a .plx file
%
% [n, ts, sv] = plx_event_ts(filename, channel)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   channel - 1-based external channel number
%             strobed channel has channel number 257  
% OUTPUT:
%   n - number of timestamps
%   ts - array of timestamps (in seconds)
%   sv - array of strobed event values (filled only if channel is 257)

[n, ts, sv] = mexPlex(3,filename, ch);