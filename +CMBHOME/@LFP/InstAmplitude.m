function signal_amplitude = InstAmplitude(signal)
% Returns instantaneous amplitude using the Hilbert transform for 'signal'
%
% signal_phase = LFP.InstAmplitude(signal);

signal_amplitude = abs(hilbert(signal));

end