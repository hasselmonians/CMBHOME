function signal_filtered = ThetaFilter(signal, Fs)
% signal_filtered = ThetaFilter(signal, Fs)
%
% Takes 'signal' and bandpasses it to theta frequencies (6 to 10 Hz)
%
% Arguments
%
% signal - arbitrary signal to be filtered
% Fs - sample frequency

% andrew 17 november 2009 adapted from tophercode

Wn_theta = [6/(Fs/2) 10/(Fs/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(3,Wn_theta);

signal_filtered = filtfilt(btheta,atheta,signal);
