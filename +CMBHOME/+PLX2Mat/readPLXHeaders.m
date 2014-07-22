function headers = readPLXHeaders(filename)
% headers = readPLXHeaders(filename)
%
% This function reads the headers from 'filename', a .plx file. It is
% important for the CMBHOME.Import.PLX function because the plexon binaries
% do not return channel numbers of event strings. This is a major bug,
% since event strings will not be propery matched up with their occurances.
%
% Ben Kraus wrote the bulk of this function and Andrew Bogaard edited the
% section to include event header information.
%
% Ben also points out: "Plexon didn't have much foresight for the number of
% units/channel, and only left enough space in the headers for 4
% units/channel. Therefore there is no way to get an accurate count of the
% number of timestamps without reading the entire file, and you should
% ignore that information stored in the header.
% 
% To put another way:
% [tscounts, wfcounts, evcounts] = plx_info(filename, fullread)
% 
% If "fullread" is false, it pulls the data straight from the header,
% which can't be relied upon to be accurate. If "fullread" is true, it
% actually reads the entire file to count how many of each type of event
% occurred. My function (readPLXHeader.m) is only designed to read the
% header, so it can't get the accurate counts."
%
% 4/1/2011

if(nargin ~= 1)
   disp('1 input arguments are required')
   return
end

if(isempty(filename))
   [fname, pathname] = uigetfile('*.plx', 'Select a plx file');
	filename = strcat(pathname, fname);
end

fid = fopen(filename, 'r');
if(fid == -1)
	disp('cannot open file');
   return
end

% read file header
headers.magic = fread(fid, 4, '*char')';
if(~strcmp(headers.magic,'PLEX'))
    error('This file is not a valid .plx file');
end

headers.version = fread(fid, 1, 'int32');
headers.comment = deblank(fread(fid, 128, '*char')');
headers.ADFrequency = fread(fid, 1, 'int32');
headers.numDSPChannels = fread(fid, 1, 'int32');
headers.numEventChannels = fread(fid, 1, 'int32');
headers.numSlowChannels = fread(fid, 1, 'int32');
headers.numPointsWave = fread(fid, 1, 'int32');
headers.numPointsPreThr = fread(fid, 1, 'int32');

YR = fread(fid, 1, 'int32');
MO = fread(fid, 1, 'int32');
DA = fread(fid, 1, 'int32');
HR = fread(fid, 1, 'int32');
MI = fread(fid, 1, 'int32');
SC = fread(fid, 1, 'int32');
headers.date = datenum([YR, MO, DA, HR, MI, SC]);

headers.fastread = fread(fid, 1, 'int32');
headers.waveformfreq = fread(fid, 1, 'int32');
headers.lasttimestamp = fread(fid, 1, 'double');

if(headers.version >= 103)
    headers.trodalness = fread(fid, 1, 'uint8');
    headers.datatrodalness = fread(fid, 1, 'uint8');
    headers.bitsperspikesample = fread(fid, 1, 'uint8');
    headers.bitsperslowsample = fread(fid, 1, 'uint8');
    headers.spikeMaxMagnitudeMV = fread(fid, 1, 'ushort');
    headers.slowMaxMagnitudeMV = fread(fid, 1, 'ushort');
else
    headers.trodalness = 1;
    headers.datatrodalness = 1;
    headers.bitsperspikesample = 12;
    headers.bitsperslowsample = 12;
    headers.spikeMaxMagnitudeMV = 3000;
    headers.slowMaxMagnitudeMV = 5000;
    fread(fid, 8, 'uint8');
end

if(headers.version >= 105)
    headers.spikePreAmpGain = fread(fid, 1, 'ushort');
else
    headers.spikePreAmyGain = 1000;
    fread(fid, 1, 'ushort');
end

if(headers.version >= 106)
    headers.acquiringsoftware = deblank(fread(fid, 18, '*char')');
    headers.processingsoftware = deblank(fread(fid, 18, '*char')');
    fread(fid, 10, 'uint8');
else
    headers.acquiringsoftware = '';
    headers.processingsoftware = '';
    fread(fid, 46, 'char');
end

headers.tscounts = fread(fid, [5, 130], 'int32');
headers.wfcounts = fread(fid, [5, 130], 'int32');
headers.evcounts = fread(fid, [1, 512], 'int32');

for ii = 1:headers.numDSPChannels;
    headers.chans(ii).name = deblank(fread(fid, 32, '*char')');
    headers.chans(ii).signame = deblank(fread(fid, 32, '*char')');
    headers.chans(ii).channel = fread(fid, 1, 'int32');
    headers.chans(ii).wfrate = fread(fid, 1, 'int32');
    headers.chans(ii).sig = fread(fid, 1, 'int32');
    headers.chans(ii).ref = fread(fid, 1, 'int32');
    headers.chans(ii).gain = fread(fid, 1, 'int32');
    headers.chans(ii).filter = fread(fid, 1, 'int32');
    headers.chans(ii).threshold = fread(fid, 1, 'int32');
    headers.chans(ii).method = fread(fid, 1, 'int32');
    headers.chans(ii).nunits = fread(fid, 1, 'int32');
    headers.chans(ii).template = fread(fid, [5, 64], 'short');
    headers.chans(ii).fit = fread(fid, 5, 'int32');
    headers.chans(ii).sortwidth = fread(fid, 1, 'int32');
    headers.chans(ii).boxes = zeros(5, 2, 4);
    headers.chans(ii).boxes(:) = fread(fid, 5*2*4, 'short');
    headers.chans(ii).sortbeg = fread(fid, 1, 'int32');
    if(headers.version>=105)
        headers.chans(ii).comment = deblank(fread(fid, 128, '*char')');
    else
        headers.chans(ii).comment = '';
        fread(fid, 128, 'char');
    end
    
    if(headers.version>=106)
        headers.chans(ii).srcid = fread(fid, 1, 'uchar');
        fread(fid, 1, 'uchar');
        headers.chans(ii).chanid = fread(fid, 1, 'ushort');
        fread(fid, 10, 'int32');
    else fread(fid, 11, 'int32');
    end
end

for ii = 1:headers.numEventChannels;
    headers.evchans(ii).name = deblank(fread(fid, 32, '*char')');
    headers.evchans(ii).channel = fread(fid, 1, 'int32');

    if(headers.version>=105)
        headers.evchans(ii).comment = deblank(fread(fid, 128, '*char')');
    else
        headers.evchans(ii).comment = '';
        fread(fid, 128, 'char');
    end
    
    if(headers.version>=106)
        headers.evchans(ii).srcid = fread(fid, 1, 'uchar');
        fread(fid, 1, 'uchar');
        headers.evchans(ii).chanid = fread(fid, 1, 'ushort');
        fread(fid, 32, 'int32');
    else fread(fid, 33, 'int32');
    end
end

fclose(fid);