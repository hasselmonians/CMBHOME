function signal_phase = InstPhase(signal)
% signal_phase = LFP.InstPhase(signal);
%
% Returns instantaneous phase using the Hilbert transform for 'signal'
%
% andrew bogaard 27 sept 2010

signal_phase = atan2(imag(hilbert(signal)), signal);

end