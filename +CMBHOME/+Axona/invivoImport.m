function session = invivoImport(varargin)
  % invivoImport(collectionFolder,overwrite,rat_folders,prefixes)
  % Takes in dacq files and creates a CMBobject out of them.
  %
  % base_path is folder containing rat folders
  % overwrite =1 to overwrite existing files, default is 0
  % set "folders" to only look in certain rats folders
  % set prefixes to only run certain prefixes of files

  % wchapman 20130927



  %% Parse inputs
  p = inputParser;
  p.addParamValue('basePath', pwd,    @(x) ischar(x));            % Look for rat folders in this folder
  p.addParamValue('cutFiles', [],      @(x) iscell(x));
  p.addParamValue('setFiles', [],     @(x) iscell(x));

  p.parse(varargin{:});
  basePath = p.Results.basePath;
  cutFiles = p.Results.cutFiles;
  setFiles = p.Results.setFiles;

  import CMBHOME.*
    
  if ~exist([basePath filesep cutFiles{1}],'file')
      error('cut file doesn''t exist');
  end
  
  %% Create Session structures:
  session = struct();
  for k = 1:length(setFiles)
      session(end+1).setFile = [basePath filesep setFiles{k}];   
      session(end).mat = [session(end).setFile(1:end-4) '.mat'];
      session(end).basePath = basePath;
      session(end).prefix = setFiles{k}(1:end-4);
      session(end).cut = cutFiles;
  end
  session(1) = [];
  session = session(:);


  %########################################################################
  %###### Now for each element in session, make the root objects###########
  %########################################################################
  startSpike = ones(size(cutFiles));
  for i = 1:length(session)
    
    %% Behavioral Data
    % Position information
    pos_file=[session(i).setFile(1:end-4) '.pos'];
    session(i).pos_file = pos_file;

    [scaleFactor, vid_ts, vid_x, vid_y, vid_headdir, d_epoch] = dacq_pos(session(i).pos_file);
    %vid_ts = vid_ts + offset;
    %offset = offset + vid_ts(end);
    vid_ts = vid_ts(:);
%     scaleFactor = 0.2326; % new default

    vid_headdir = wrapTo180(vid_headdir);

    %create root structure
    root = CMBHOME.Session('name', (session(i).mat), ...
      'b_ts', vid_ts, ...
      'b_x', vid_x, ...
      'b_y', vid_y, ...
      'raw_pos', 1, ...
      'raw_headdir', 1, ...
      'b_headdir', vid_headdir, ...
      'date_created', now, ...
      'epoch', [-inf inf], ...
      'fs_video', 1/mean(diff(vid_ts)), ...
      'spatial_scale',scaleFactor); % put it all together
    
    %% spike structure
    for k = 1:length(cutFiles)
        if isempty(strfind(cutFiles{1},'clu'))
            inds = strfind(cutFiles{k},'_');
            tetNum = str2num(cutFiles{k}(inds(end)+1:end-4));
            dotnums{k} = [basePath filesep session(i).prefix '.' num2str(tetNum)];
        else
            inds = strfind(cutFiles{k},'.');
            tetNum = str2num(cutFiles{k}(inds(end)+1:end));
            dotnums{k} = [basePath filesep session(i).prefix '.' num2str(tetNum)];
        end
    end
    
    [spike, spikeWaveforms, ~, startSpike] = dacq_tetrode(basePath, cutFiles, dotnums, startSpike);

    root.spike = spike;
    
    % Change the waveform format into that expected by Visualize2 instead
    % of Ehren's visualization method
    clear waveform
    waveform(1,4).mean = [];
    waveform(1,4).std = [];
    
    
    for m = 1:size(spikeWaveforms,2)
        for p = 1:size(spikeWaveforms,1)
            fid = fopen(session(i).setFile);
            query = ['gain_ch_' num2str(p) ' '];
            s = textscan(fid, '%s', 'delimiter', '\n');
            idx = find(arrayfun(@(x) ~isempty(find(strfind(x{1}, query))),s{1}));
            gain = str2num(s{1}{idx}(length(query):end));
            gain = (3000000/gain)/256;
            fclose(fid);
            
            if ~isempty(spikeWaveforms(p,m).mean)
                for q = 1:4
                    if q==1
                        waveform(end+1,q).mean = spikeWaveforms(p,m).mean(:,q)*gain;
                    else
                        waveform(end,q).mean = spikeWaveforms(p,m).mean(:,q)*gain;
                    end
                    waveform(end,q).std = spikeWaveforms(p,m).std(:,q)*gain;
                end
            end
        end
    end
    waveform(1,:) = [];
    
    
    %{
    for k = 1:4
        waveform(:,k).mean = arrayfun(@(x) x.mean(:,k), spikeWaveforms,'UniformOutput',0);
        waveform(:,k).std = arrayfun(@(x) x.std(:,k), spikeWaveforms,'UniformOutput',0);
    end
    waveform=waveform(:);
    
    for k = 1:size(waveform)
        for m = 1:4
            w(k,m).mean = waveform(k).mean{m}(:)';
            w(k,m).std = waveform(k).std{m}(:)';
        end
    end
    %}
    root.user_def.waveform=waveform; 
    
    %align
    root = root.AlignSpike2Session;
     
    %% LFP:
    try
        path_lfp = dacq_eeg(basePath, session(i).prefix);
        root.path_lfp = path_lfp(:);
    end
    
    root.user_def.dacqInfo = dacq_set(session(i).setFile);
    root.user_def.importInfo = session(i);
    
    %% Save it
    save(session(i).mat,'root');
    
  end
  
end
