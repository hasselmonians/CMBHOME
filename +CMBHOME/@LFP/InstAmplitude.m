function signal_amplitude = InstAmplitude(signal)
% signal_phase = LFP.InstAmplitude(signal);
%
% Returns instantaneous amplitude using the Hilbert transform for 'signal'
%
% eric zilli 3 nov 2010

signal_amplitude = abs(hilbert(signal));

end