function [adfreq, n, ad] = plx_ad_span_v(filename, ch, startCount,endCount)
import CMBHOME.PLX2Mat.*
% plx_ad_span_v(filename, channel): Read a span of a/d data from a .plx file
%
% [adfreq, n, ad] = plx_ad_span_v(filename, ch,startCount,endCount)
%
% INPUT:
%   filename - if empty string, will use File Open dialog
%   startCount - index of first sample to fetch
%   endCount - index of last sample to fetch
%   channel - 0 - based channel number
%
% OUTPUT:
%   adfreq - digitization frequency for this channel
%   n - total number of data points 
%   ad - array of a/d values converted to mV

[adfreq, n, ad] = mexPlex(18,filename, ch, startCount, endCount);
