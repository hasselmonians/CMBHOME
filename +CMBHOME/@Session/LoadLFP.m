function self = LoadLFP(self, ind, varargin)
% Loads LFP data file into root.b_lfp property
%
% (1) Loads root.path_lfp(root.active_lfp)
% (2) Loads all root.path_lfp(indices) to current object
%
% So far, .mat files with lfp, lfp_ts fields can be imported. .ncs
% files can be imported. other methods
% can be added
%
% (1) root = root.LoadLFP
% (2) root = root.LoadLFP(indices)

%% Setup
import CMBHOME.Utils.* %need path2lfp

if ~exist('ind', 'var')
  if isempty(self.active_lfp), error('root.LoadLFP: Please specify root.active_lfp');
  else ind = self.active_lfp; end
elseif isempty(ind)
  warning('CMBH:error', 'root.LoadLFP: Indeces in root = root.LoadLFP(indeces) was empty, so nothing was done. Do not include argument indeces to use root.active_lfp.');
  return
end


if ~iscell(self.path_lfp), error('root.LoadLFP: Please specify LFP files via root.AddLFP'); end
if ~isa(self.b_lfp, 'CMBHOME.LFP'), self.b_lfp = CMBHOME.LFP; end


p = inputParser;
p.addParamValue('ifScale',  1,           @(x) (x==1) || (x==0));           % scale eeg data?
p.addParamValue('ifHD', 1,               @(x) (x==1) || (x==0));           % use egf instead of eeg if available?
p.addParamValue('downsample', 0,         @(x) isnumeric(x));               % [0] downsample to default of 250 Hz if set to 1, or to specific sampling rate if set to anything other than 1

p.parse(varargin{:});
ifScale = p.Results.ifScale;
ifHD = p.Results.ifHD;
downsample = p.Results.downsample;


for i = 1:length(ind)
  f_name = self.path_lfp{ind(i)};
  %%  If is cell array then make absolute path
  if iscell(f_name),
    if exist(CMBHOME.Utils.array2path(f_name,dataPath),'file')
      f_name = CMBHOME.Utils.array2path(f_name,dataPath); %f_name = f_name(2:end);
    elseif exist(fullfile(dataPath,CMBHOME.Utils.array2path(f_name,dataPath)),'file')
      f_name = fullfile(dataPath,array2path(self.path_lfp{ind(i)}));
    else
      error('couldnt find specified LFP file: %s',CMBHOME.Utils.array2path(f_name));
    end
  end
  
  if exist(f_name, 'file')
    ext = lower(f_name(find(f_name=='.', 1, 'last'):end));
    ext = ext(1:4);
    
    switch ext
      
      %-----------------------------------------------------------------%
      %-----------------------------.EEG--------------------------------%
      %-----------------------------------------------------------------%
      case '.eeg'
        
        % if HD requested and possible, update pathLFP and recall LoadLFP
        HDf_name = [f_name(1:find(f_name=='.', 1, 'last')+1),'gf',f_name(find(f_name=='.', 1, 'last')+4:end)];
        if ifHD && exist(HDf_name,'file')
          eegTxtInd = strfind(self.path_lfp{ind(i)}{end},'.eeg');
          self.path_lfp{ind(i)}{end}(eegTxtInd+1:eegTxtInd+3) = 'egf';
          if isempty(varargin)
            self = self.LoadLFP(ind(i));
          else
            self = self.LoadLFP(ind(i),varargin{:});
          end
        else
          %otherwise load lowdef eeg
          fid = fopen(f_name,'r');
          if (fid == -1)
            error('Could not open file %s',f_name);
          end
          for ii = 1:8
            textstring = fgetl(fid);
          end
          fs = sscanf(textstring,'%*s %f');
          for ii = 1:3
            textstring = fgetl(fid);
          end
          nosamples = sscanf(textstring,'%*s %u');
          fseek(fid,10,0);
          signal = fread(fid,nosamples,'int8');
          ts = (0:length(signal)-1)'/fs;
          fclose(fid);
          
          channel_name = f_name;
          
          % Check to see if we need to reverse the polarity
          try %#ok<*TRYNC>
            method = self.user_def.dacqInfo.eeg(ind).mode;
            if strcmp(method,'-sig')
              signal = -1*signal;
              disp('Reversed sign of LFP signal')
            end
          end
          
          % seal the deal
          self.b_lfp(ind(i)) = CMBHOME.LFP(signal, ts, fs, channel_name);
          
          % downsample if requested
          if downsample & self.b_lfp(ind(i)).fs > downsample
            ts_old = self.b_lfp(ind(i)).ts;
            ts_new = linspace(ts_old(1),ts_old(end),round(length(ts_old)*(downsample/self.b_lfp(ind(i)).fs)))';
            self.b_lfp(ind(i)).signal = interp1(ts_old,self.b_lfp(ind(i)).signal,ts_new);
            self.b_lfp(ind(i)).fs = downsample;
            self.b_lfp(ind(i)).ts = ts_new;
            fprintf('Downsampled to %i Hz\n',self.b_lfp(ind(i)).fs);
          end
          
          % Scale the eeg (to mV) if desired
          if ifScale == 1
            scaleFactor = 11.7188; % @ gain of 1000x min to max range is -1.5 to 1.5 mV, loaded EEG ranges from -128 to 127 in integer steps with bit precision, thus this value was calculated as 3000uV / 256steps to give 11.7188 uV/step.
            gain = self.user_def.dacqInfo.eeg(ind(i)).gain;
            self.b_lfp(ind(i)).signal = scaleFactor * (self.b_lfp(ind(i)).signal/gain);
          end
          
        end
        
        %-----------------------------------------------------------------%
        %-----------------------------.MAT--------------------------------%
        %-----------------------------------------------------------------%
        
      case '.mat'
        
        load(f_name)
        
        channel_name = ind(i);
        
        if exist('lfp', 'var')
          signal = lfp;
        else
          warning('CMBH:error', 'Vector with signal not found.')
        end
        
        if exist('lfp_ts', 'var')
          ts = lfp_ts;
          fs = mean(diff(lfp_ts))^-1;
        else
          warning('CMBH:error', 'Vector with signal timestamps not found. Cannot assign root.fs')
        end
        
        self.b_lfp(ind(i)) = CMBHOME.LFP(signal, ts, fs, channel_name);
        
        
        %-----------------------------------------------------------------%
        %-----------------------------.PLX--------------------------------%
        %-----------------------------------------------------------------%
        
      case '.plx'
        % not written yet
        warning('CMBH:error', 'Not coded up yet');
        
        %-----------------------------------------------------------------%
        %-----------------------------.EGF--------------------------------%
        %-----------------------------------------------------------------%
      case '.egf'
        
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
        
        % Check to see if we need to reverse the polarity
        try
          method = self.user_def.dacqInfo.eeg(ind).mode;
          if strcmp(method,'-sig')
            signal = -1*signal;
            disp('Reversed sign of LFP signal')
          end
        end
        
        % seal the deal
        self.b_lfp(ind(i)) = CMBHOME.LFP(signal, ts, fs, channel_name);
        
        % downsample if requested
        if downsample & self.b_lfp(ind(i)).fs > downsample
          ts_old = self.b_lfp(ind(i)).ts;
          ts_new = linspace(ts_old(1),ts_old(end),round(length(ts_old)*(downsample/self.b_lfp(ind(i)).fs)))';
          self.b_lfp(ind(i)).signal = interp1(ts_old,self.b_lfp(ind(i)).signal,ts_new);
          self.b_lfp(ind(i)).fs = downsample;
          self.b_lfp(ind(i)).ts = ts_new;
          fprintf('Downsampled to %i Hz\n',self.b_lfp(ind(i)).fs);
        end
        
        % Scale the eeg (to mV) if desired
        if ifScale == 1
          scaleFactor = 0.0458; % 3000uV range at 1000x gain conversion calculated for 16 bit precision used in 4800Hz recordings.
          gain = self.user_def.dacqInfo.eeg(ind(i)).gain;
          self.b_lfp(ind(i)).signal = scaleFactor * (self.b_lfp(ind(i)).signal/gain);
        end
        
        
        %-----------------------------------------------------------------%
        %-----------------------------.NCS--------------------------------%
        %-----------------------------------------------------------------%
      case '.ncs'
        % 'n_samples' below is a vector with an integer indicating
        % how many samples exist at a certain index that
        % you can pass at the fifth parameter to NlxMatCSC
        import CMBHOME.NL2Mat.*
        
        %import_segment_length = 2^14; % each segment has 2^12 samples
        
        [fs, ~] = Nlx2MatCSC(f_name, [0 0 1 1 0], 0, 1); % check the sampling frequency, and get the number of samples
        %     1. Timestamps
        %     2. Channel Numbers
        %     3. Sample Frequency
        %     4. Number of Valid Samples
        %     5. Samples
        
        if length(unique(fs))~=1, warning('CMBH:notify', '%s', ['Sampling rate in ' fname ' is not consistent. Check timestamps.']); end
        
        decimation_factor = round(fs(1)/400); % we want to downsample to around 400 Hz
        
        decimation_factor = 2^(nextpow2(decimation_factor)-1); % simplifies our downsampling between time and signal
        
        fs_new = fs(1)/decimation_factor;
        
        warning('CMBH:notify', '%s', ['Downsampled from ' num2str(fs(1)) 'Hz to ' num2str(fs_new) 'Hz.']);
        
        [tmp_ts, samples] = Nlx2MatCSC(f_name, [1 0 0 0 1], 0, 1);
        
        samples = samples(:);
        signal = samples(1:decimation_factor:length(samples));
       
        l_seg = floor(length(signal)/length(tmp_ts));
        
        ts = zeros( length(tmp_ts) * l_seg , 1 );
        
        tmp_ts(end+1) = tmp_ts(end)+mean(diff(tmp_ts)); %#ok<AGROW>
        
        for ii = 1:length(tmp_ts)-1
          
          dt = (tmp_ts(ii+1)-tmp_ts(ii))/l_seg * (0:l_seg);
          dt = dt(1:end-1);
          ts( (ii-1) * l_seg + 1 : ii * l_seg) = tmp_ts(ii) + dt;
          
        end
        
        ts = ts*10^-6; % convert from microseconds to seconds
        
        % reshape the lfp signal
        ts = ts(:);
        channel_name = self.path_lfp{ind(i)};
        fs = fs_new;
        
        self.b_lfp(ind(i)) = CMBHOME.LFP(signal, ts, fs, channel_name);
        
        
    end
    
    
  else
    warning('CMBH:error', ['It appears file ' f_name ' in root.path_lfp does not exist, and was skipped. Run root.AddLFP to find LFP file(s) again.']);
  end
end

self = self.AlignSpike2LFP;

end
