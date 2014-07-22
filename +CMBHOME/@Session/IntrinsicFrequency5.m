function [f_peak, thetaMod, fullResults] = IntrinsicFrequency5(self, cel, print_to_screen, speed_thresh, varargin)
% [f_peak, thetaMod, rVal, fullResults] = IntrinsicFrequency3(self, cel, [print_to_screen], [speed_thresh], [varargin])
%
% ARGUMENTS
%
%   cel                     [tetrode index, cell index] in Session object root
%
%   [print_to_screen]       optional - 0 or 1 (0) plots analysis       
%
%   [speed_thresh]          optional - [min velocity, max velocity] ([0 5]
%                           cm/sec). If [-1 -1], no speed threshold
%
%   [params]                optional input params: 
%                             'printtext',bool
%                             'subtract_mean',bool
%                             'acorr',{cor,lag}
%                             'freqs',[lowBound highBound]
%
% RETURNS
%
%   f_peak                  scalar, frequency peak within theta range
%                           (empty if not theta modulated cell)
%
%   thetaMod                scalar, computed as peak(peak1 - trough1) / peak1 
%                           where peak1 is the amplitude of the first peak 
%                           that is not centered at t=0 and trough1 is the
%                           amplitude of the trough between 0 and peak1.
%                           Note: this is only computed if the correlation
%                           between the cosine wave with frequency f_peak
%                           is significantly correlated with the
%                           autocorrelation, otherwise it defaults to nan.
%
%   fullResults             struct, all values computed that may be of interest 
%                          
%
% If root.IntrinsicFrequency is called with multiple epochs, the analysis
% is performed on each epoch separately, and f_peak and power_ratio are
% vectors of length "number of epochs".
%
% THE ALGORITHM
%
% Calculate unbiased xcorr of spike trains that occur during epochs of 
% running longer than .5 seconds. For each frequency over the range of 4 to 
% 12 Hz, compute the dot product between the xcorr and a cosine wave
% starting at the first zero crossing and continuing to the third zero
% crossing (comprising a full cycle) and then standardize the result by the
% dot product between the xcorr and the cosine wave between the 2nd and 3rd
% zero crossings (comprising the peak). f_peak is taken to the frequency of the wave with the
% highest resulting value.  

% modified from intrinsicFrequency3 on 120706 by eln

% applicable?
%If the correlation is above 0.162
% (the value at which p<0.05 with df=100), the depth of modulation is
% computed to be the difference between the amplitude of the first peak and
% the first trough divided by the amplitude of the first peak.  
%
% eln 20130811


default_speed_thresh = 5; %cm/sec

p = inputParser;

p.addRequired('self');
p.addRequired('cel', @isnumeric);
p.addRequired('print_to_screen');
p.addOptional('speed_thresh', [0 default_speed_thresh], @isnumeric);
p.addParamValue('printtext', 1, @(c) numel(c)==1);
p.addParamValue('subtract_mean', 0, @(c) c==1 || c==0)
p.addParamValue('acorr', [], @(c) iscell(c) && all(size(c)==[1,2]))
p.addParamValue('freqs', 4:.01:12, @isnumeric);
p.addParamValue('version', 'firstCycle', @ischar);
p.addParamValue('downsample', 0, @isnumeric)

p.parse(self, cel, print_to_screen, speed_thresh, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
print_to_screen = p.Results.print_to_screen;
speed_thresh = p.Results.speed_thresh;
printtext = p.Results.printtext;
subtract_mean = p.Results.subtract_mean;
acorr = p.Results.acorr;
freqs = p.Results.freqs;
version = p.Results.version;
downsample = p.Results.downsample;

import CMBHOME.Utils.*

f_peak = []; % initialize vars
theta_modulated = [];

pow_rat_thresh = 1.5;

if ~exist('speed_thresh', 'var')% if no velocity threshold set, make upper bound=default, lower=0
    
    speed_thresh = [0 default_speed_thresh];
    
elseif length(speed_thresh) ~= 2
    
    disp('Improper speed thresholding');
    return
    
elseif diff(speed_thresh)<0
    
    disp('Improper speed thresholding');
    return
    
end

if ~exist('print_to_screen', 'var'), print_to_screen = 0; end

if isempty(acorr)
	orig_epoch = self.epoch;
else
	orig_epoch = [0 inf];
end

f_peak = zeros(size(orig_epoch, 1), 1);

thetaMod = zeros(size(orig_epoch, 1), 1);

acN = zeros(size(orig_epoch, 1), 1);

for epoch_ind = 1:size(orig_epoch, 1)
    
	if isempty(acorr)
    self.epoch = orig_epoch(epoch_ind,:);

    max_lag_t = .6; % 600 ms long x-corr 

    t_bin = .005; % 1 ms bin width   
    
    if downsample == 0
      
      %[cor, lag] = self.AutoCorr(cel, 'speed_thresh', speed_thresh, 'theta_skip', 0, 't_bin', t_bin, 'supress_plot', 1, 'max_lag', max_lag_t,'unbiased',1);
      
      [cor, lag, acN(epoch_ind)] = downsample4AC(self.spk.ts, t_bin, max_lag_t, downsample);
      
    else
      
      if isnan(downsample), [f_peak, thetaMod, fullResults] = deal(nan); return; end
      
      [cor, lag, acN(epoch_ind)] = downsample4AC(self.spk.ts, t_bin, max_lag_t, downsample);
      
    end
  
        
	else
        
        % Use precomputed autocorrelation
        cor = acorr{1}; 
        lag = acorr{2};
				
				if any(lag<0)
					cor = cor(lag>0);
					lag = lag(lag>0);
				end
				
				max_lag_t = max(lag);
				
				t_bin = mean(diff(lag));
				
	end
	
	if subtract_mean
    
      cor = cor - mean(cor);
      
      pow_rat_thresh = 0; % this is unknown for now
    
	end
	
  switch(version)
		case 'firstCycle'
			
			n = length(cor);
			freqs = freqs(1.5.*(1./freqs) <= max(lag) );
			
			dPhi = (freqs' .* (2*pi)) ./ t_bin^-1;
			phases = cumsum(repmat(dPhi,1,n) .* ones(length(freqs),n),2);
			wave = cos(phases);
			
			% build full cycle mask
			fullCycleMask = phases>(1*pi) & phases<=(3*pi);
            
      % cut out only portion of interest
      fullCycleWave = wave; fullCycleWave(~fullCycleMask) = nan;

			% build peak mask
			peakMask = phases>(3*pi/2) & phases<=(5*pi/2);
			
% 			% fix area under the positive portion of each curve to be 1
% 			fullCycleWave = fullCycleWave ./ repmat(nansum(fullCycleWave.*peakMask,2),1,n);
			
			% remove oversampled frequencies
			widths = sum(peakMask,2);
      [~,uniqInds] = unique(widths);
			freqs = freqs(uniqInds);
      fullCycleWave = fullCycleWave(uniqInds,:);

			R = corr(cor', fullCycleWave','rows','pairwise');

			% this method given relative dip
% 			% 			fullCycleSpkRate = dot(fullCycleMask.*wave,repmat(cor',length(freqs),1),2);
% 			for f = 1:length(freqs)
% 				x = fullCycleWave(f,:);
% 				y = cor;
% 				y(isnan(x)) = [];
% 				x(isnan(x)) = [];
% 				fullCycleSpkRate(f) = dot(x,y);
% 			end
% 			peakSpkRate = dot(peakMask.*wave,repmat(cor',length(freqs),1),2);
% 			R = fullCycleSpkRate./peakSpkRate;
% 			plot(lag,cor)
% 			hold on
% 			plot(1./freqs,(R/2 +0.5) .*max(cor),'-r')
% 			hold off
% 			pause
		case 'wavelet'
			
			ac = reshape(cor,1,[]) - mean(cor);
			ac = [fliplr(ac),ac];
			R = nan(1,length(freqs));
			for fi = 1:length(freqs)
				width = 6;
				
				dt = mean(diff(lag));
				sf = freqs(fi)/width; %f/width;
				st = 1/(2*pi*sf);
				
				t=-3.5*st:dt:3.5*st;
				m = morlet(freqs(fi),t,width);
				
				ac_mid = floor(length(ac)/2);
				m_half = floor(length(m)/2);
				if mod(length(m),2)
					tmpAC = ac(ac_mid-m_half:ac_mid+m_half);
				else
					tmpAC = ac(ac_mid-m_half+1:ac_mid+m_half);
				end
				R(fi) = abs(dot(tmpAC,m));
			end
			R = R ./ sum(R);
			
	end
	
	
	[~,imax] = extrema(R); % find local maxima and see if any lie within theta
	ftmp = freqs(imax); ftmp = reshape(ftmp,[],1);
	stmp = R(imax); stmp = reshape(stmp,[],1);
	stmp = stmp(ftmp > 6 & ftmp < 9);
	ftmp = ftmp(ftmp > 6 & ftmp < 9);
	
	[~,ind] = max(stmp); % intrinsic frequency
	if ~isempty(ind)  % if theta peak was found
		f_peak(epoch_ind) = ftmp(ind);
		thetaMod(epoch_ind) = stmp(ind);
	else                                                % if theta peak was not found
		f_peak(epoch_ind) = nan;
		thetaMod(epoch_ind) = nan;
	end
	
	
	fullResults(epoch_ind).lag = lag;
	fullResults(epoch_ind).cor = cor;
	fullResults(epoch_ind).freqs = freqs;
	fullResults(epoch_ind).f_i = ind;
	fullResults(epoch_ind).f_peak = f_peak;
	fullResults(epoch_ind).f_thetaMod = thetaMod;
	fullResults(epoch_ind).modFactors = R;
	
	
    if print_to_screen
			figure
			subplot(1,2,1)
			plot(freqs,R); hold on; A = axis;
			plot([f_peak(epoch_ind) f_peak(epoch_ind)], A(3:4),'k','linewidth',2);
			xlabel('Frequency (Hz)')
			ylabel('Theta Modulation')
      ys = ylim;
      
			if thetaMod>0
			txt{1} =  ['Intrinsic Frequency = ', num2str(f_peak(epoch_ind)), ' Hz'];
			txt{2} =  ['Modulation depth = ', num2str(thetaMod(epoch_ind))];
				text(diff(xlim)/2, ys(2)*.90,txt, 'FontSize', 8);
      elseif thetaMod(epoch_ind) < pow_rat_thresh && printtext
        text(diff(xlim)/2, ys(2)*.90, 'Not theta-modulated.', 'FontSize', 8);
      end

      hold off
      
			subplot(1,2,2)
			f_wave = (wave(ind,:)-min(wave(ind,:)))/diff(minmax(wave(ind,:)));
			f_wave = f_wave.*diff(minmax(cor)) + min(cor);
			bar(lag,cor);
			hold on;
			plot(lag,f_wave,'r','linewidth',2);
    end
    
end

end
