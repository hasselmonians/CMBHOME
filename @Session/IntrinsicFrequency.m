function [f_peak, power_ratio] = IntrinsicFrequency(self, cel, print_to_screen, speed_thresh, varargin)
% (1) F = root.IntrinsicFrequency(cel);
% (2) F = root.IntrinsicFrequency(cel, print_to_screen)
% (3) [F, power_ratio] = root.IntrinsicFrequency(self, cel, print_to_screen, speed_thresh)
%
% ARGUMENTS
%
%   cel                     [tetrode index, cell index] in Session object root
%   print_to_screen         0 or 1 (0) plots analysis       
%   speed_thresh            [min velocity, max velocity] ([0 5] cm/sec). If [-1 -1], no
%                           speed threshold
%   params                  (see below)
%
% RETURNS
%
%   F                       scalar, frequency peak within theta range
%                           (empty if not theta modulated cell)
%   power_ratio             scalar, average power within 1 Hz of the peak
%                           in the theta band over total signal power
%
% If root.IntrinsicFrequency is called with multiple epochs, the analysis
% is performed on each epoch separately, and f_peak and power_ratio are
% vectors of length "number of epochs".
%
% THE ALGORITHM
%
% See Jeewajee et al, 2008 Hippocampus.
%
% Calculate unbiased xcorr of spike trains that 
% occur during epochs of running longer than .5 seconds, take fft 
% of the xcorr, show the power from 0 - 15, find peak in the 7-11Hz band. 
% If the mean power within 1Hz of the peak is 50% greater than the mean 
% spectral power (calculated as the mean of the square of the signal),
% then the cell is theta modulated
%
% andrew bogaard 12 july 2010

default_speed_thresh = 5; %cm/sec

p = inputParser;

p.addRequired('self');
p.addRequired('cel', @isnumeric);
p.addRequired('print_to_screen');
p.addOptional('speed_thresh', [0 default_speed_thresh], @isnumeric);
p.addParamValue('figure_handle', '', @(c) numel(c)==1);
p.addParamValue('std_smooth_kernel', 3, @isnumeric);
p.addParamValue('printtext', 1, @(c) numel(c)==1);
p.addParamValue('chronux_spectrum', 1, @(c) c==1 || c==0)
p.addParamValue('subtract_mean', 0, @(c) c==1 || c==0)

p.parse(self, cel, print_to_screen, speed_thresh, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
print_to_screen = p.Results.print_to_screen;
speed_thresh = p.Results.speed_thresh;
figure_handle = p.Results.figure_handle;
std_smooth_kernel = p.Results.std_smooth_kernel;
printtext = p.Results.printtext;
chronux_spectrum = p.Results.chronux_spectrum;
subtract_mean = p.Results.subtract_mean;

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

power_ratio = zeros(size(orig_epoch, 1), 1);

for epoch_ind = 1:size(orig_epoch, 1)
    
    self.epoch = orig_epoch(epoch_ind,:);

    max_lag_t = .6; 

    t_bin = .002; % 2 ms bin width   
    
    [cor, lag] = self.AutoCorr(cel, 'speed_thresh', speed_thresh, 'theta_skip', 0, 't_bin', t_bin, 'supress_plot', 1);

    if subtract_mean
    
      cor = cor - mean(cor);
      
      pow_rat_thresh = 0; % this is unknown for now
    
    end
    
    if chronux_spectrum, [S, f] = getChronuxSpectrum(cor, t_bin); % use chronux if user requests
    else [S, f] = getMatlabSpectrum(cor, t_bin); % use matlab by default  
    end

    [xmax,imax] = extrema(S); % find local maxima and see if any lie within theta

    ftmp = f(imax);
    stmp = S(imax);

    [a_peak,ind] = max(stmp(ftmp > 7 & ftmp < 11)); % intrinsic frequency

    ftmp = ftmp(ftmp > 7 & ftmp < 11);

    if ~isempty(ind)                                    % if theta peak was found
        f_peak(epoch_ind) = ftmp(ind)';     
        a_peak_av = mean(S(f>f_peak(epoch_ind)-1 & f<f_peak(epoch_ind)+1));
        peak = S(f==f_peak(epoch_ind));
    else                                                % if theta peak was not found
        f_peak(epoch_ind) = 0;
        peak = 0;
        a_peak_av = 0;
    end

    power_ratio(epoch_ind) = a_peak_av / mean(S);

    if print_to_screen
        
        line(f(f>=0 & f<=15), 10^9*S(f>=0 & f<=15), 'Color', 'k', 'LineWidth', 2), hold on;

        ymax = max([10^9*peak*1.1 10^9*1.5*mean(S)/4]); 

        if ymax>0, ylim([0 ymax]); end
        
        xlim([0 15]);

        ylabel('Power (x 10^-9)')
        xlabel('Frequency (hz)');

        ys = ylim;

        line([f_peak(epoch_ind) f_peak(epoch_ind)], ys, 'Color', 'r', 'linewidth', 1.5), hold on;

        if power_ratio(epoch_ind) >= pow_rat_thresh && printtext
            text(f_peak(epoch_ind)+15*.02, ys(2)-diff(ys)*.04, ['Intrinsic Frequency = ' num2str(f_peak(epoch_ind)) ' Hz'], 'FontSize', 8);
        elseif power_ratio(epoch_ind) < pow_rat_thresh && printtext
            text(8, ys(2)-diff(ys)*.04, 'Not theta-modulated.', 'FontSize', 8);
        end

        hold off
        
    end
    
end

end

function [S, f] = getChronuxSpectrum(cor, t_bin)
% CHRONUX spectrum of correlation

    import CMBHOME.Utils.*

    AddChronuxPackage; % add chronux functionality

    params.tapers = [1.25, 1];
    params.fpass = [];
    params.Fs = 1/t_bin;
    params.pad = 6;  

    [S,f]=mtspectrumc(cor,params); % compute spectrum of cor

    % stdf = .2/mean(diff(f)); % elements in a .2 hz std
    % 
    % kernel = pdf('Normal', -3*stdf:3*stdf, 0, stdf); % smooth with gaussian kernel

    % S = conv(S, kernel, 'same'); % plot(f, S) 

    % S = 10*log10(S); %convert to decibels

end

function [S, f] = getMatlabSpectrum(cor, t_bin)
% MATLAB

    % padd array with zeros to 2^16 elements
    if length(cor)<2^16, cor = cat(1, zeros(floor((2^16-length(cor))/2),1), cor); end
    if length(cor)<2^16, cor = cat(1, cor, zeros(2^16-length(cor),1)); end

    S = (abs(fft(cor)).^2) / length(cor);       %Compute the spectrum.
    S = S(1:length(cor)/2+1);                   %Ignore negative frequencies.
                                                %Define the time resolution.
    f0 = 1/t_bin;                               %Determine the sampling frequency.
    df = 1/(t_bin*length(cor));                 %Determine the frequency resolution.
    fNQ = f0 / 2;                               %Determine the Nyquist frequency.
    f = (0:df:fNQ);                             %Construct frequency axis.

    % stdf = .2/mean(diff(f)); % elements in a .2 hz std
    % kernel = pdf('Normal', -3*stdf:3*stdf, 0, stdf); % smooth with gaussian kernel
    % 
    % %S = conv(S, kernel, 'same'); % plot(f, S) 
    % 
    % S = 10*log10(S/max(S));            %Convert to decibels.

end
