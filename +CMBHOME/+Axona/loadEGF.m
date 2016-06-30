function signal = loadEGF(f_name, ifScale)
    fid = fopen(f_name,'r');
    if (fid == -1)
      error('Could not open file %s',f_name);
    end
    for ii = 1:8
      textstring = fgetl(fid);
    end
    fs = sscanf(textstring,'%*s %f');
    for ii = 1:2
      textstring = fgetl(fid);
    end
    nosamples = sscanf(textstring,'%*s %u');
    fseek(fid,10,0);
    signal = fread(fid,nosamples,'int16');
    ts = (0:length(signal)-1)'/fs;
    fclose(fid);

    channel_name = f_name;

    % Scale the eeg (to mV) if desired
    if ~exist('ifScale','var')
        ifScale = 0;
    end
    
    if ifScale == 1
      scaleFactor = 0.0458; % 3000uV range at 1000x gain conversion calculated for 16 bit precision used in 4800Hz recordings.
      gain = 1;
      signal = scaleFactor * (signal/gain);
    end
    
end
        