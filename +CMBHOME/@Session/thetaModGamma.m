function [modindex, thetarange, gammarange, powPhsDists, bincenters, partial] = thetaModGamma(self,varargin)
% Computes the phase modulation of gamma by theta 
% using the method described by Tort et al (2010) J Neurophys.
% In short, it computes the
% entropy in the histogram of powers over phases.
%
%
% INPUT ARGS
%  root - cmb object with active_lfp set and lfp loaded
%  filtType - default = [], otherwise 'theta' or 'vel' to process data only
%             with high theta power or fast running speeds.
%  filtParams - default = [], lower and upper bounds to use in the
%               filtering of the data as mentioned above
%  shuffles - default 0, use if you want to compute null dist of modulation
%  thetaRange - default = [4:0.25:12] frequencies to look at phases of 
%  gammaRange - defualt = [30:1:120] frequencies to look at the power of
%  ifPlot - Default = 0. If 1, plots in new figure. If anything else, then
%           assumes input is an axis to plot onto
%
%
% OUTPUT ARGS
%  modindex - matrix of size [nG x nT] where nT is the number of theta frequencies 
%             analyzed and nG is the number of gamma frequencies analysed
%  thetarange - same as inputed, usable as labels for the matrix when plotting
%  gammarange - same as inputed, usable as labels for the matrix when plotting
%  powPhsDists - normalized power histograms over phases for theta band that 
%                showed the max modulation. has shape [nG x nT x nPh] where
%                nT is the number of theta frequencies analyzed and nG is
%                the number of gamma frequencies analysed and nPh is the
%                number of phase bins used.
%  bincenters - phase bin labels
%
% [modindex, thetarange, gammarange, powPhsDists, bincenters] = thetaModGamma(root,filtType,filtParams,shuffles,thetarange,gammarange)


    p = inputParser;
    p.addParamValue('filtType', []);
    p.addParamValue('filtParams',   []);
    p.addParamValue('shuffles', 0);
    p.addParamValue('thetarange',   4:0.25:12);
    p.addParamValue('gammarange',   30:1:120);
    p.addParamValue('ifPlot',   0);

    p.parse(varargin{:});
    
    thetarange = p.Results.thetarange;
    gammarange = p.Results.gammarange;
    filtType = p.Results.filtType;
    filtParams = p.Results.filtParams;
    shuffles = p.Results.shuffles;
    ifPlot = p.Results.ifPlot;
   
  epochs = self.epoch;
  self.epoch = [0 inf];
    
  % initialize outputs
  bincenters = nan(36,1);

  if ~exist('partial','var')||isempty(partial)
      % get gamma amplitudes
      [~,gammaamps] = multiphasevec(gammarange,self.lfp.signal,self.lfp.fs,8);
      
      % get theta filtered signals and phase
      thetaangles = extractThetaPhase(self,'wavelet',thetarange);
      partial.gammaamps = gammaamps;
      partial.thetaangles = thetaangles;
  else
      gammaamps = partial.gammaamps;
      thetaangles = partial.thetaangles;
  end
  
  % select out high quality data (no saturations and no flat lines + apply requested filtering
  [gammaamps,thetaangles] = cleanData(self,gammaamps,thetaangles,thetarange,filtType,filtParams,epochs);

  % double check that there are at least 200 cycles of theta (7Hz), the min
  % mentioned by Tort et al (2010).
  if size(thetaangles,2)<200*(self.lfp.fs/8)
    warning('There may be too few theta cycles!');
  end

  if ~shuffles
    % compute the values of interest
    [modindex, powPhsDists, bincenters] = computeModIndex(gammaamps,thetaangles,thetarange);  
  else
    % compute significance thresholds for each frequency pairing
    [modindex, powPhsDists, bincenters] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange);
  end
  
  if ifPlot
     if ifPlot==1
        figure;
        ax = gca;
     else
        ax = ifPlot; 
     end
     
     axes(ax);
     imagesc(thetarange,gammarange,modindex);
     axis xy
     
     title(['thetaModGamma  (maxmod = ' num2str(max(modindex(:)))])
     xlabel('Theta Frequency'), ylabel('Gamma Frequency')
  end
  
end



%% Cleandata
function [gammaamps,thetaangles] = cleanData(root,gammaamps,thetaangles,thetarange,filtType,filtParams,epochs)

  % compute specific filters if needed
  if exist('filtType','var') && ~isempty(filtType)

    % find good amps
    switch(filtType)
      case 'theta'
        thetaamps = abs(hilbert(buttfilt(root.lfp.signal,minmax(thetarange),root.lfp.fs,'bandpass',4)));
        excInds = thetaamps < max(thetaamps)*filtParams(1) | thetaamps > max(thetaamps)*filtParams(2);
      case 'vel'
        hdvel = interp1(root.ts,root.vel,root.lfp.ts);
        excInds = hdvel < filtParams(1) | hdvel > filtParams(2);
      otherwise
        error('Unknown filtType');
    end
    excInds = reshape(excInds,1,[]);
  else
    excInds = zeros(1,size(gammaamps,2));
  end

  %%  drop saturated or flatlined data
  
  %wchapman: Dropped these because they are specific to the recording system
  %maxSig = 126-1e-5; % threshold for railed out at top 
  %minSig = -127+1e-5; % threshold for railed out at bottom
  
  maxSig = max(root.lfp.signal);
  minSig = min(root.lfp.signal);
  t0 = [0-1e-5 0+1e-5]; % thresholds for no dV/dt
  
  badInds = [...
    CMBHOME.Utils.ThresholdBandDetect(diff(root.lfp.signal),t0(1),t0(2),1,5);... % flatline data
    CMBHOME.Utils.OverThresholdDetect(root.lfp.signal,maxSig,1,5);... % railed out at top
    CMBHOME.Utils.UnderThresholdDetect(root.lfp.signal,minSig,1,5)... % railed out at bottom
    ];

  if 0;%~isempty(badInds)
    OE = root.epoch;
    % add one second buffer around each bit of bad data to remove potential
    % edge effects on the pahse estimates
    %badInds = root.lfp.ts(min(max([badInds(:,1)-root.lfp.fs badInds(:,2)+root.lfp.fs],1),length(root.lfp.ts)));
    if isvector(badInds), badInds = reshape(badInds,1,[]); end
    % many windows probably overlap now, merge them
    badInds = CMBHOME.Utils.MergeEpochs2(badInds,1);
    % reapply whatever window the original epoching was
    badInds = CMBHOME.Utils.IntersectEpochs2(badInds,root.epoch);
    % define epochs as only the bad data with the surrounding buffers
    if ~isempty(badInds)
        root.epoch = badInds;
        % build an 'inds' vector in lfp object to find inds of the bad data
        root.b_lfp(root.active_lfp).b_myvar = 1:length(root.b_lfp(root.active_lfp).signal);
        % initialize vector for data to remove
        badInds = zeros(1,size(gammaamps,2));
        % use inds into bad data epochs to set the boolean flags
        badInds(CMBHOME.Utils.ContinuizeEpochs(root.lfp.myvar)) = 1;
        % restore epoch to normal so we can still refer to root.lfp
        root.epoch = OE;
    else
        badInds = zeros(1,size(gammaamps,2));
    end
  else
    badInds = zeros(1,size(gammaamps,2));
  end

  %% apply epoch
  root.epoch = [0 inf];
  outepoch = root.lfp.ts;
  root.epoch = epochs;
  if size(epochs,1)>1
      try
        outepoch(ismember(outepoch,cell2mat(root.lfp.ts')))=NaN;
      catch
         outepoch(ismember(outepoch,cell2mat(root.lfp.ts)))=NaN; 
      end
  else
    outepoch(ismember(outepoch,root.lfp.ts'))=NaN;
  end
  outepoch = ~isnan(outepoch);
  root.epoch = [0 inf];

  %%
  
  % perform filter
  gammaamps(:,excInds|badInds|outepoch') = [];
  thetaangles(:,excInds|badInds|outepoch') = [];

end


%% 
function [modindex, powPhsDists, bincenters] = computeModIndex(gammaamps,thetaangles,thetarange)
  % bin gamma amplitudes by theta phases
  meanamps = nan(36,size(gammaamps,1),size(thetaangles,1));
  bincenters = nan(36,size(thetaangles,1));
  for i=1:length(thetarange)
    [meanamps(:,:,i),bincenters]=Findmeanamps(thetaangles(i,:)',gammaamps');
  end
  
  % get back to the way things used to be
  meanamps = permute(meanamps,[2 3 1]);

  % normalize the area to one
  meanampsNorm = meanamps ./ repmat(sum(meanamps,3),[1,1,size(meanamps,3)]);  

  % compute the mod indexes
  shannonent = sum(meanampsNorm .* log10(meanampsNorm),3);
  modindex = (log10(36)+shannonent) ./ log10(36);

  % ELN edit 120512 - The following code will change powPhsDists to be the
  % power distrobution over theta phases for the theta band with the
  % maximum modulation index.
  %%% compute mean phase alignments over window of high modulation
  powPhsDists = meanamps;
end

%%
function [modindex_out, powPhsDists_out, bincenters] = computeShuffledModIndex(shuffles,gammaamps,thetaangles,thetarange)
  fprintf('Shuffling (%i): ',shuffles);
  modindex = nan(size(gammaamps,1),size(thetaangles,1),shuffles);
  powPhsDists = nan(size(gammaamps,1),36,shuffles);

  for i = 1:shuffles
    % keep user feel confident that all is well
    fprintf('%i ',i);
    
    % find random offset to shift all theta angles by
    offsetOkay = false;
    while ~offsetOkay
      offset = round(rand*size(thetaangles,2));
      offsetOkay = offset>1 && offset<size(thetaangles,2);
    end

    % shift theta angles
    ta = [thetaangles(:,offset+1:end), thetaangles(:,1:offset)];
    % compute mod indices for this shift
    [modindex(:,:,i), powPhsDists(:,:,i), bincenters] = computeModIndex(gammaamps,ta,thetarange);
  end

  % find distro values
  modindex_out(:,:,1) = mean(modindex,3);
  modindex_out(:,:,2) = std(modindex,[],3);
  powPhsDists_out(:,:,1) = mean(powPhsDists,3);
  powPhsDists_out(:,:,2) = std(powPhsDists,[],3);
  fprintf('\n');
end

%%
function [meanamps bincenters]=Findmeanamps(slowphase,fastamp)
  %this function finds the normalized mean amplitude of the fast wave for each phase of the
  %slow wave
  %it recieves the phase of the slow wave (in radians from 0-2pi) and the
  %corresponding amplitude of the fast wave (in mV)

  %it returns a matrix with the normalized mean amplitude (row 2) corresponding to the
  %phases(row 1-- max of corresponding phase), and the middle of the phase
  %for a bar graph(row 3)

  
  
  %%make into column vectors if not
  if size(slowphase,1)==1
    slowphase=slowphase';
  end

  % double check data matches
  if size(fastamp,1)~=size(slowphase,1)
    error('fastamp must have time in rows');
  end

  % define phase bins
  sorting = linspace(-pi, pi, 37);
  
  % perform sorting
  meanamps = nan(length(sorting)-1,size(fastamp,2));
  for bin = 1:length(sorting)-1
    currInds = slowphase > sorting(bin) & slowphase <= sorting(bin+1);
    meanamps(bin,:) = mean(fastamp(currInds,:),1);
  end
  
  % add in bin labsls
  bincenters = mean([sorting(1:end-1);sorting(2:end)])'+pi;

end

%%
function [phase,pow]=multiphasevec(f,S,Fs,width)
% FUNCTION [phase,pow]=multiphasevec(f,S,Fs,width)
%
% Returns the phase and power as a function of time for a range of
% frequencies (f).
%
% Simply calls phasevec in a loop.
%
% INPUT ARGS:
%   f = [2 4 8];   % Frequencies of interest
%   S = dat;       % Signal to process
%   Fs = 256;      % Sampling frequency
%   width = 6;     % Width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   phase- Phase data [freqs,time]
%   power- Power data [freqs,time]
%

pow = nan(length(f),length(S)); 
phase = nan(length(f),length(S));

for a=1:length(f)
     [phase(a,:),pow(a,:)]=phasevec(f(a),S,Fs,width);
end

end

%%
function thetaPhs = extractThetaPhase(root,method,band)

if isempty(root.lfp)
	error('Please load lfp with active_lfp before calling');
end

switch(lower(method))
	
	case 'waveform'
		if exist('band','var') || ~isempty(band)
			warning('band not used when computing phase by waveform')
		end
		
		%  method described by 
		broadBand = [1 60];
		filtSig = buttfilt(root.lfp.signal,broadBand,root.lfp.fs,'bandpass',4);
		
		% define window detection params
		min_sep = root.lfp.fs / 20; % 50ms separation between peak and next trough
		min_length = min_sep; % same for the separation between trough and next peak
		
		% Find troughs and peaks 
		% as onsets and offsets of signal with positive slope
		trphsAndPeaks = CMBHOME.Utils.OverThresholdDetect(sign(diff(filtSig)),0,min_sep,min_length);
		
		% Find zero crossings as sign changes in signal
		upAndDownZeroXings = CMBHOME.Utils.OverThresholdDetect(sign(filtSig),0,min_sep,min_length);

		waveLandmarks = sort([trphsAndPeaks, upAndDownZeroXings]);
		
		error('This code is not complete')
		% this method was to be based off of the methods described by Belluscio
		% et al 2012, however they had a large number of channels
		% simultaneously recorded from which they computed the median signal
		% and then computed the wavefrom from this.  This likely dramatically
		% reduces the number of local minima and maxima found and thereby makes
		% this a practical endeavor.  To make this work on a signal channel of
		% data, the signal would likely have to be filtered more thereby
		% reducing the value of using this method.
		
	case 'wavelet'
		% returns a vector of values, spanning theta range, per time point 
		if exist('band','var') & ~isempty(band)
			thetarange = band;
		else
			thetarange = 4:0.1:12;
		end
		[thetaPhs,~] = multiphasevec(thetarange,root.lfp.signal,root.lfp.fs,16);

		
	case 'hilbert'
		% returns a single phase estimate by time point using bandpass filtering
		if exist('band','var') & ~isempty(band)
			thetarange = band;
		else
			thetarange = [4 12];
		end
		thetaPhs = angle(hilbert(buttfilt(root.lfp.signal,thetarange,root.lfp.fs,'bandpass',4)));
	
	otherwise
		error('Unknown method')
end
end

%%
function y = morlet(f,t,width)
% function y = morlet(f,t,width)
% 
% Morlet's wavelet for frequency f and time t. 
% The wavelet will be normalized so the total energy is 1.
% width defines the ``width'' of the wavelet. 
% A value >= 5 is suggested.
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)
%
% See also: PHASEGRAM, PHASEVEC, WAVEGRAM, ENERGY 
%
% Ole Jensen, August 1998 

sf = f/width;
st = 1/(2*pi*sf);
A = 1/sqrt(st*sqrt(pi));
y = A*exp(-t.^2/(2*st^2)).*exp(i*2*pi*f.*t);
end

%%
function [y, pow] = phasevec(f,s,Fs,width)
% FUNCTION [y,pow] = phasevec(f,s,Fs,width)
%
% Return the phase as a function of time for frequency f. 
% The phase is calculated using Morlet's wavelets. 
%
% INPUT ARGS:
%   f = 8;     % Frequency of interest
%   s = dat;   % The eeg signal to evaluate
%   Fs = 256;  % Sample rate
%   width = 6; % width of Morlet wavelet (>= 5 suggested).
%
% OUTPUT ARGS:
%   y- Phase as a function of time
%   pow- Power as a function of time
%
% Ref: Tallon-Baudry et al., J. Neurosci. 15, 722-734 (1997)

dt = 1/Fs;
sf = f/width;
st = 1/(2*pi*sf);

t=-3.5*st:dt:3.5*st;
m = morlet(f,t,width);

y = conv(s,m);

pow = abs(y).^2;
pow = pow(ceil(length(m)/2):length(pow)-floor(length(m)/2));

l = find(abs(y) == 0); 
y(l) = 1;

% normalizes phase estimates to length one
y = y./abs(y);
y(l) = 0;
   
y = angle( y(ceil(length(m)/2):length(y)-floor(length(m)/2)) );
end
