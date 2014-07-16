function signal_filtered = DeltaFilter(signal, Fs)

%% Takes 'signal' and bandpasses it for delta frequencies (2 to 4 Hz)

% Arguments

% signal - arbitrary signal to be filtered
% Fs - sample frequency (specific to your recording rig) ex. Mark's
% Neuralynx system has a sampling period of .000528, or 1893.9393... Hz

% andrew 17 november 2008 adapted from tophercode

Wn_theta = [2/(Fs/2) 4/(Fs/2)]; % normalized by the nyquist frequency

[btheta,atheta] = butter(3,Wn_theta);

signal_filtered = filtfilt(btheta,atheta,signal);
