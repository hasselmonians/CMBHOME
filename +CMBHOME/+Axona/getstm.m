function [stimts] = getstm(datafile)

% Read STM file
%
% Usage: [stimts] = getstm(datafile);
%
%
%
% Jim Donnett <donnett@axona.com>
% Axona Ltd
% 148 Green Lane, St. Albans, Herts, AL3 6EU, U.K.
% http://www.axona.com
%
% This script makes hard-coded assumptions about the .stm file format, namely that
% the timebase is on header line 8, and the number of samples on line 11, followed
% by data_start on line 12.
% 
% Copyright (C) 2005  Axona Ltd
% All Rights Reserved
%
% This M-file is released under Q Public License v 1.0,
% a copy of which should accompany this file.

fid = fopen(datafile,'r');
if (fid == -1)
   error(sprintf('Could not open file %s',filename));
end
for i = 1:8
   textstring = fgetl(fid);
end
Fs = sscanf(textstring,'%*s %f');
for i = 1:3
   textstring = fgetl(fid);
end
nosamples = sscanf(textstring,'%*s %u');
fseek(fid,10,0);
stimts = fread(fid,nosamples,'uint32',0,'b');
fclose(fid);
% divide timestamps by timebase to convert into seconds
stimts = stimts/Fs;
