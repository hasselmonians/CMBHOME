function waveStruct = ReadDotNum(fname)
% waves = ReadDotNum('C:\Users\Holger\Desktop\Test data Ehren\CRE45\140513_1432_pedestalCircleTrack.2')

packetLength = 216; %Given in DacqUSBFileFormats.PDF

%open raw event file and error check
fid = fopen(fname,'r','ieee-be'); %file formats PDF specifies big-endian

if (fid == -1)
    warning('Cannot open file %s, will place hold',file);
else
    
    for k = 1:14, fgetl(fid); end %14 full lines of metadata
    
    headerSize = ftell(fid) + 10; %skip "data_start" then immediately begin binary
    fseek(fid,0,1); % find end of file
    Npackets = ((ftell(fid)-headerSize)/packetLength); % determine total number of packets in file
    
    fseek(fid,headerSize,-1); %rewind the file to end of header
    
    % exit smoothly if failed to read in an integer number of packets
    if ~isinteger(Npackets), Npackets = floor(Npackets); end;
    
    %For each packet, read the timestamp and 50 8bit samples for all 4 channels
    time_stamps = nan(1,Npackets);
    wave1 = int8(nan(50,Npackets));
    wave2 = int8(nan(50,Npackets));
    wave3 = int8(nan(50,Npackets));
    wave4 = int8(nan(50,Npackets));
    
    for j = 1:Npackets
        time_stamps(j)= fread(fid,1,'int32')/96000; %divide by time base
        wave1(:,j)=fread(fid,50,'int8');
        
        [~] = fread(fid,1,'int32')/96000;
        wave2(:,j)=fread(fid,50,'int8');
        
        [~] = fread(fid,1,'int32')/96000;
        wave3(:,j)=fread(fid,50,'int8');
        
        [~] = fread(fid,1,'int32')/96000;
        wave4(:,j)=fread(fid,50,'int8');
    end
    
    waveStruct.ts = time_stamps;
    waveStruct.t1 = double(wave1);
    waveStruct.t2 = double(wave2);
    waveStruct.t3 = double(wave3);
    waveStruct.t4 = double(wave4);
    
end