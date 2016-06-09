function [f_peak, thetaMod, rVal, acN, fullResults] = IntrinsicFrequency3(self, cel, print_to_screen, speed_thresh, varargin)
% Intrinsic frequency of a unit
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
%                             'printtext'
%                             'subtract_mean'
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
%   rVal                    scalar, the correlation btw the cosine wave and the xcorr
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
% See Sharp et al 2008 Hippocampus.
%
% Calculate unbiased xcorr of spike trains that occur during epochs of
% running longer than .5 seconds. Compute the correlation between the xcorr
% and each of a set of cosine waves with frequencies between 4-12Hz in
% 0.01Hz steps.  f_peak is taken to the frequency of the wave with the
% highest correlation with the xcorr.  If the correlation is above 0.162
% (the value at which p<0.05 with df=100), the depth of modulation is
% computed to be the difference between the amplitude of the first peak and
% the first trough divided by the amplitude of the first peak.
%
% [f_peak, thetaMod, rVal, fullResults] = IntrinsicFrequency3(self, cel, [print_to_screen], [speed_thresh], [varargin])

default_speed_thresh = 5; %cm/sec

p = inputParser;

p.addRequired('self');
p.addRequired('cel', @isnumeric);
p.addRequired('print_to_screen');
p.addOptional('speed_thresh', [0 default_speed_thresh], @isnumeric);
p.addParamValue('printtext', 1, @(c) numel(c)==1);
p.addParamValue('subtract_mean', 0, @(c) c==1 | c==0)
p.addParamValue('downsample', 0, @isnumeric)

p.parse(self, cel, print_to_screen, speed_thresh, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
print_to_screen = p.Results.print_to_screen;
speed_thresh = p.Results.speed_thresh;
printtext = p.Results.printtext;
subtract_mean = p.Results.subtract_mean;
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

orig_epoch = self.epoch;

f_peak = zeros(size(orig_epoch, 1), 1);

rVal = zeros(size(orig_epoch, 1), 1);

thetaMod = zeros(size(orig_epoch, 1), 1);

acN = zeros(size(orig_epoch, 1), 1);

max_lag_t = .5; % 500 ms long x-corr

for epoch_ind = 1:size(orig_epoch, 1)
  
  self.epoch = orig_epoch(epoch_ind,:);
  
  t_bin = .005; % 1 ms bin width
  
  if downsample == 0
    
    %[cor, lag] = self.AutoCorr(cel, 'speed_thresh', speed_thresh, 'theta_skip', 0, 't_bin', t_bin, 'supress_plot', 1, 'max_lag', max_lag_t,'unbiased',1);
    
    [cor, lag, acN(epoch_ind)] = downsample4AC(self.cel_ts, t_bin, max_lag_t, downsample);
    
  else
    
    if isnan(downsample), [f_peak, thetaMod, rVal, fullResults] = deal(nan); return; end
    
    [cor, lag, acN(epoch_ind)] = downsample4AC(self.cel_ts, t_bin, max_lag_t, downsample);
    
  end
  
  if subtract_mean
    
    cor = cor - mean(cor);
    
    pow_rat_thresh = 0; % this is unknown for now
    
  end
  
  freqs = 4:.01:12;
  wave = nan(length(freqs),length(cor));
  for f = 1:length(freqs)
    wave(f,:) = cos(linspace(0,freqs(f)*max_lag_t*2*pi,length(cor)));
  end
  R = corr(cor(:),wave');
  
  [rVal(epoch_ind),ind] = max(R);
  
  f_peak(epoch_ind) = freqs(ind);
  
  Rthresh = 0.164; % as used in Sharp & Koester (2008)
  
  % compute 1st-peak minus trough
  smoothWdw = 25; %ms
  smthCor = smooth(cor, ceil(smoothWdw/(t_bin*1000)) );
  [~,iMax,~,iMin] = extrema(wave(ind,:));
  iMax = sort(iMax);
  iMin = sort(iMin);
  % replace first max with real max (cos always 1 at 0)
  [~,iMax(1)] = max(smthCor(iMax(1):iMin(1)));
  % find window to compute mean over
  wndw = ceil((25/2)/(t_bin*1000)); %25 ms window
  % find trough
  [~, trough_i] = min(smthCor(iMax(1):iMax(2)));
  trough_i = trough_i + iMax(1) - 1;
  trough_mag = mean(smthCor(max(trough_i-wndw,1):trough_i+wndw));
  % find peak
  [~,peak_i] = max(smthCor(iMin(1):iMin(2)));
  peak_i = peak_i + iMin(1) - 1;
  peak_mag = mean(smthCor(peak_i-wndw:peak_i+wndw));
  
  thetaMod(epoch_ind) = (peak_mag - trough_mag) / peak_mag;
  
  fullResults(epoch_ind).t = lag;
  fullResults(epoch_ind).cor = cor;
  fullResults(epoch_ind).freqs = freqs;
  fullResults(epoch_ind).R = R;
  fullResults(epoch_ind).f_i = ind;
  fullResults(epoch_ind).f_peak = freqs(ind);
  fullResults(epoch_ind).f_R = R(ind);
  fullResults(epoch_ind).smoothedCor = smthCor;
  fullResults(epoch_ind).trough_mag = trough_mag;
  fullResults(epoch_ind).trough_i = trough_i;
  fullResults(epoch_ind).peak_mag = peak_mag;
  fullResults(epoch_ind).peak_i = peak_i;
  fullResults(epoch_ind).thetaMod = thetaMod(epoch_ind);
  fullResults(epoch_ind).rhythmic = rVal(epoch_ind) > Rthresh;
  
  if print_to_screen
    
    t = fullResults(epoch_ind).t;
    plot(t, fullResults(epoch_ind).cor, 'Color', 'k', 'LineWidth', 2); hold on;
    
    plot(t, wave(ind,:)*std(cor)+mean(cor),'r','LineWidth',2);
    
    if ~isnan(thetaMod)
      A = axis;
      plot(t(peak_i), cor(peak_i),'.g','markersize',20);
      plot(t(trough_i), cor(trough_i),'.g','markersize',20);
      plot([A(1) A(2)], [cor(peak_i) cor(peak_i)], '--g');
      plot([A(1) A(2)], [cor(trough_i) cor(trough_i)], '--g');
    end
    
    
    ys = ylim;
    
    if ~isnan(thetaMod) && printtext
      txt{1} =  ['Intrinsic Frequency = ', num2str(f_peak(epoch_ind)), ' Hz'];
      txt{2} =  ['rho = ',num2str(rVal(epoch_ind))];
      txt{3} =  ['Modulation depth = ', num2str(thetaMod)];
      text(diff(xlim)/2, ys(2)*.94,txt, 'FontSize', 8);
    elseif power_ratio(epoch_ind) < pow_rat_thresh && printtext
      text(diff(xlim)/2, ys(2)*.94, 'Not theta-modulated.', 'FontSize', 8);
    end
    
    hold off
    
  end
  
end

end
