function ImportGUISaddoris(ifExternal)

    if ~exist('ifExternal','var')
        ifExternal = 1;
    end

    addpath(genpath('~/github/CMBHOME'))
    import CMBHOME.*
    
    %% Graphically Cue for the input files
    [vf, vp, ~] = uigetfile('*', 'VideoFile');
    videoFile = [vp vf];
    [f, p, ~] = uigetfile('*', 'SpikeFile');
    spikeFile = [p f];
    [f, p, ~] = uigetfile('*', 'LFPFile');
    lfpFile = [p f];
        
    
    %% New Video Load
    % Note, this uses python. Ask Bill for "ImageReader" environment for
    % list of modules to install.
    
    % For help on calling python from matlab, see "'Call User-Defined
    % Python Module" in the matlab 'doc' interface
    
    % To call from MATLAB:
    if ifExternal == 0
        cmd = ['/anaconda/envs/ImageReader/bin/python ReadBehavior.py "' videoFile '" 0'];
        unix(cmd)
        
    else
        % Easier and faster to call from Terminal!
        exec = '/anaconda/envs/ImageReader/bin/python';
        cmd = ['\n Run this: \n\n' exec ' ReadBehavior.py "' videoFile '" 1' '\n\n'];
        
        fprintf(cmd)
        
        input('Press any key once Python program complete: ...')
    end
    
    % After running, import to MATLAB
    posFile = [videoFile(1:end-3) 'csv'];    
    delimiter = ' ';    
    formatSpec = '%f%f%[^\n\r]';
    fileID = fopen(posFile,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    
    y = dataArray{:, 1};
    x = dataArray{:, 2};
    y = max(y) - y;

    vid_ts = (cumsum(ones(size(x)))-1)/30;
    
    %% Initialize the structure with just behavioral variables 
    root = CMBHOME.Session('b_x',x,'b_y',y,...
                            'b_headdir', zeros(size(x)),... % ... for now
                            'fs_video', 30,...
                            'b_ts', vid_ts);

    %% Import the Spikes
    plx = CMBHOME.PLX.readPLXFileC(spikeFile,'all');
    for i = 1:numel(plx.SpikeChannels)
        for k = 0:max(plx.SpikeChannels(i).Units) 
            
                inds = plx.SpikeChannels(i).Units==k;
                ts = plx.SpikeChannels(i).Timestamps(inds);
                ts = double(ts) / 40000;
                waves = plx.SpikeChannels(i).Waves(:,inds);
                waves = double(waves);
                
                if ~isempty(ts)
                    if k == 0
                        stashind = max(plx.SpikeChannels(i).Units)+1;
                    else
                        stashind = k;
                    end

                    waveform(i,stashind).mean = mean((waves),2)';
                    waveform(i,stashind).std = std((waves),[],2)';


                    Spike(i,stashind) = CMBHOME.Spike('ts',ts,'vid_ts',vid_ts);
                end


        end
    end

    root.spike = Spike;
    root.user_def.waveform = waveform;

    %% new import LFP?
    for chan = 1:16
        [adfreq, n, ts, fn, ad] = plx_ad_v(lfpFile, chan+16); % Odd offset, but seems consistent in plx files
        ts = (cumsum(ones(size(ad)))-1)/adfreq;
        fs = adfreq;
        
         % Downsampling
        decimation_factor = round(fs(1)/1000); % we want to downsample to around 400 Hz
        decimation_factor = 2^(nextpow2(decimation_factor)-1); % simplifies our downsampling between time and signal
        fs_new = fs(1)/decimation_factor;
        
        sig = ad(1:decimation_factor:length(ad));
        ts = ts(1:decimation_factor:length(ts));

        % Import
        LFP(chan) = CMBHOME.LFP(sig,ts,fs_new,[]);
        
    end
    root.b_lfp = LFP;
    
    %% Fix it up
    root = root.FixTime;
    root = root.FixPos;
    
    CMBH = [videoFile(1:end-3) 'mat'];
    save(CMBH,'root');
    
end