function [ts, data] = dacq_inp(fname)

%% Header
fid = fopen(fname,'r','b');
fseek(fid,0,-1);

for i = 1:15
    header{i} = fgetl(fid);
end
timebase = str2num(header{10}(10:end-3));
numPackets = str2num(header{15}(17:end));

% "datastart"
for i = 1:10
    fread(fid,1,'uint8=>char');
end
dataStart = ftell(fid);
%% Each packet
ts = fread(fid,numPackets,'int32',3); 
fseek(fid,dataStart,-1); fread(fid,1,'int32');
type = fread(fid,numPackets,'uint8=>char',6);
fseek(fid,dataStart,-1); fread(fid,1,'int32'); fread(fid,1,'uint8');

for i = 1:16
    fseek(fid,dataStart,-1); fread(fid,1,'int32'); fread(fid,1,'uint8');
    fread(fid,i-1,'ubit1');
    data(:,i) = fread(fid,numPackets,'ubit1',55);
end

%% Process for output
bads = find(type~='I');

ts(bads) = [];
data(bads,:) = [];

ts = ts / timebase;

end