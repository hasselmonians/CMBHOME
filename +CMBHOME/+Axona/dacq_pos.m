function [scaleFactor vid_ts vid_x vid_y headdir d_epoch] = dacq_pos(fname)
%[root.scaleFactor, root.ts, root.x, root.y, root.headdir, root.epoch] =
%dacq_pos(dataPath,prefix)
%
% Finds video file with the prefix "prefix" in folder "dataPath", and
% returns relevant spatial information for a CMBobject

% Bill 2012.06.20

%file information
%file = strcat(prefix,'.pos');
%dataFile = fullfile(dataPath,file);
import CMBHOME.*
packetLength = 20; % pos packet is only 20 bytes 
% format=
% 1- 4 4bytes = frame counter, 
% 5- 6 2bytes = bigLED X pos
% 7- 8 2bytes = bigLED Y pos
% 9-10 2bytes = smallLED X pos
%11-12 2bytes = smallLED Y pos
%13-14 2bytes = number_of_pixels_in_big_spot
%15-16 2bytes = number_of_pixels_in_little_spot
%17-18 2bytes = total_tracked_pixels
%19-20 2bytes = empty

%open file and error check
%fid = fopen(dataFile,'r','ieee-be');
fid = fopen(fname,'r','ieee-be');
if (fid == -1)
  warning('Cannot open file. Please check',file);
else
    %% End time and sampling rate
    for i = 1:4, fgetl(fid); end
    endTime=fgetl(fid);
    for i = 6:17, fgetl(fid); end
    sampling_rate=fgetl(fid);
    
    [~, endTime] = strtok(endTime);
    [endTime, ~] = strtok(endTime);
    [~, sampling_rate] = strtok(sampling_rate);
    [sampling_rate, ~] = strtok(sampling_rate);
    endTime=str2double(endTime);
    sampling_rate=1/(str2double(sampling_rate));
    d_epoch=[0 endTime];
    
    % Get the spatial resolution
    for i = 19:25; fgetl(fid);end
    [~, scaleFactor] = strtok(fgetl(fid));
    scaleFactor = scaleFactor(2:end);
    scaleFactor = 100/(str2num(scaleFactor));
    
    
    %% skip past rest of header
    for i = 27; fgetl(fid);end
    headerSize = ftell(fid) + 10; 
    fseek(fid,0,1); % find end of file
    Npackets = floor((ftell(fid)-headerSize)/packetLength); % determine total number of packets in file

    dataStart=headerSize; %+ 4; % to set data start to begin after frame counter
    fseek(fid,dataStart,-1); %rewind the file to end of header
    
    % allocate 
    bx = nan(1,Npackets);
    by = nan(1,Npackets);
    lx = nan(1,Npackets);
    ly = nan(1,Npackets);
    vid_ts = nan(1,Npackets);
    
    % Loop through for each packet and return as loc_data
    for i = 1:Npackets
        vid_ts(i) = fread(fid,1,'int32');
        bx(i) = fread(fid,1,'int16');
        by(i) = fread(fid,1,'int16');
        lx(i) = fread(fid,1,'int16');
        ly(i) = fread(fid,1,'int16');
        fseek(fid,8,0);
    end
    
    fclose(fid);

    [vid_ts, i] = unique(vid_ts);
    bx = bx(i);
    by = by(i);
    lx = lx(i);
    ly = ly(i);

    %% extract X,Y for each element
    
    % 1023 stands for missing datapoint
    bx(bx==1023) = nan;
    by(by==1023) = nan;
    lx(lx==1023) = nan;
    ly(ly==1023) = nan;

    % This would be necessary to match thetint output
    %{
    bx = bx/scaleFactor;
    lx = lx/scaleFactor;
    by = by/scaleFactor;
    ly = ly/scaleFactor;
    %}
    
    ly = -ly + 1023;
    by = -by + 1023;

    %Take average of big & little dots for each time point
    vid_x = nanmean([bx;lx],1);
	vid_y = nanmean([by;ly],1);
    %vid_ts = vid_ts * sampling_rate;
    vid_ts = (1:length(vid_x)) * sampling_rate;
    vid_x = nanInterp(vid_x);
    vid_y = nanInterp(vid_y);

    
    
    %% Head direction
    
    %Find head direction from 4 quadrant arctangent
    ratio = unitsratio('deg','rad');
    headdir = ratio .* atan2( by-ly , bx-lx);
    headdir = mod(headdir,360)';		
		
    headdir = fixHeaddir(headdir);
 
end

end


