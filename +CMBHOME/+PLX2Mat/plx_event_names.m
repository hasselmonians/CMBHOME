function [n,names] = plx_event_names(filename)
import CMBHOME.PLX2Mat.*
% plx_event_names(filename): Read name for each event type from a .plx file
%
% [n,names] = plx_event_names(filename)
%
% INPUT:
%   filename - if empty string, will use File Open dialog

% OUTPUT:
%   names - array of event name strings
%   n - number of channels

[n,names] = mexPlex(16,filename);