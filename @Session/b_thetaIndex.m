function [ind,peak,cor,lag] = b_thetaIndex(self,cel,varargin)

%% Constructed from the methods from
% Grid cells without theta oscillations in the entorhinal cortex of bats
% Michael M. Yartsev,	 Menno P. Witter	 & Nachum Ulanovsky
% Nature 479, 103�107 (03 November 2011) doi:10.1038/nature10583

% Jason Climer

p = inputParser;

p.addRequired('self');
p.addRequired('cel', @isnumeric);
p.addParamValue('t_bin',0.01); % Bin size of all the temporal autocorrelograms was 10 ms  

p.parse(self, cel, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
t_bin = p.Results.t_bin;

max_lag = 0.5; % The autocorrelogram was computed for � 500 ms lags

% Acor - taken from intrinsic frequency 2
if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')
    max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
end
[cor, lag] = self.spk_xcorr(cel, max_lag, t_bin, 1, 'prob');


% The peak of the autocorrelogram at zero lag was equalized to the maximal 
% value not including the zero lag peak


% PERSONAL CORRESTPONANCE
% Regarding the ylimits of the autocorrelograms - we followed the 
% convention used in some of Moshe Abeles's papers, where Abeles used 
% autocorrelograms of spike trains (Abeles, BTW, introduced the notion of 
% autocorrelations to Neuroscience in the 1970's together with George 
% Gerstein - and we all know Abeles well here in Israel... :)  ) -- and, 
% following Abeles, we normalized the min and max of the autocorrelograms 
% to be between -1 and +1

cor = cor/max(cor(lag~=0))*2-1;

% The power spectrum of the temporal autocorrelograms was 
% assessed by computing the fast Fourier transform (FFT) of the 
% autocorrelogram, and calculating the square of the FFT magnitude
S = abs(fft(cor)).^2;
df = 1/range(lag); fNQ = 1/mode(diff(lag))/2;
f = 0:df:fNQ;
S = S(1:length(f));

% The power spectrum was smoothed with a 2-Hz rectangular window
S = smooth(S,2/df);

% and the peak value in the 5-11 Hz band was identified
peak = f(f>=5&f<=11);
[~,i] = max(S(f>=5&f<=11));
peak = peak(i);

% A neuron was defined as theta)modulated if the mean power within 1)Hz of 
% each side of the peak in the 5�11 Hz frequency range was at least 5 times
% greater than the mean spectral power between 0 Hz and 50 Hz
ind = mean(S(abs(f-peak)<1))/mean(S(f<50));

end
