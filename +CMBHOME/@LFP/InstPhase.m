function signal_phase = InstPhase(signal)
% Returns instantaneous phase using the Hilbert transform for 'signal'
%
% signal_phase = LFP.InstPhase(signal);

signal_phase = atan2(imag(hilbert(signal)), signal);

end