function signal_filtered = BandpassFilter(signal, Fs, Fpass)
% signal_filtered = BandpassFilter(signal, Fs, Fpass)
%
% Takes 'signal' and bandpasses it to Fpass frequencies
%
% Arguments
%
% signal - arbitrary signal to be filtered
% Fs - sample frequency
% Fpass - 1x2 vector of F_low and F_high indicating passband

% andrew 17 november 2009

if exist('Fpass', 'var')
    if diff(Fpass)<=0
        help CMBHOME.LFP.BandpassFilter
        error('Fpass must be in format [f_low, f_high]');
    end
else
    help CMBHOME.LFP.BandpassFilter
    error('You must pass argument Fpass. See help above.');
end        

Wn_theta = [Fpass(1)/(Fs/2) Fpass(2)/(Fs/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(3,Wn_theta);

signal_filtered = filtfilt(btheta,atheta,signal);