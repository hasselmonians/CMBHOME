function ripple_events = RippleEvents(self, ind)
% Detects onset, center, and offset of ripple events that occur in LFP
% signal. Returns indices in root.b_lfp(ind).signal
%
% Filters to ripple band (80-250Hz), finds ripple envelope signal (defined 
% by smoothed magnitude of hilbert transform) and detects threshold
% crossings. As per Nyugen paper in 2009:
%
% "For the amplitude based method, the detection signal was the ripple
% envelope, which was computed as the absolute value of the Hilbert
% transformer on the ripple-band signal (80-250 Hz), and then smoothed 
% with a 50 ms Gaussian window. [T]he upper detection threshold was set to 
% mean(detect signal) + 3 � std(detect signal), and the lower threshold 
% was set to be mean(detect signal) + 1.5 � std(detect signal). 
%
% When the ripple envelope exceeded the 
% upper threshold, this was flagged a ripple event. The extent of the ripple 
% event was determined by the first crossing of the lower threshold. For both 
% detection methods, events of less than 30 ms were discarded, and events 
% that overlapped in time with non SWS periods (determined by high 
% theta/delta power ratio) were also discarded."
%
% arguments:
% ind - index in root.b_lfp(ind).signal
%
% returns:
% ripple_events - Nx3 array of ripple event indices: [onsets, offsets, centers; ...]
%
% andrew 7 december 2009
% andrew 18 may 2010 updated for CMBHOME

import CMBHOME.*

gaussian_window = gausswin(round(.050*self.b_lfp(ind).fs)); % 50 ms gaussian window

min_duration = .020; % seconds, minimum ripple duration

ripple_envelope = sqrt( self.BandpassLFP(ind, 'Passband', 'ripple') ); % ripple envelope

ripple_envelope = conv(ripple_envelope, gaussian_window, 'same'); % convolve with gaussian

lthresh = ripple_envelope > 1.5*std(ripple_envelope)+mean(ripple_envelope); % upper threshold crossings

onsets = find(diff(lthresh)==1);
offsets = find(diff(lthresh(onsets(1):end))==-1)+onsets(1);

e = length(offsets);

l_inds = cat(2, onsets(1:e), offsets); % lower threshold crossings indeces

r_i = 1; % initialize variables
ripple_events = [];

for i = 1 : size(l_inds, 1)
    
    [maxr, j] = max(ripple_envelope(l_inds(i,1):l_inds(i,2)));
    
    if maxr > 3*std(ripple_envelope)+mean(ripple_envelope) % if the peak crosses upper threshold
        
        ripple_events(r_i, 1) = l_inds(i,1); % ripple onset
        ripple_events(r_i, 2) = l_inds(i,2); % ripple offset
        ripple_events(r_i, 3) = l_inds(i,1)+j-1; % ripple center
        
        r_i = r_i+1;
        
    end
    
end

ripple_events(ripple_events(:,2)-ripple_events(:,1)<min_duration*self.b_lfp(ind).fs,:) = []; % remove ripples with lower crossing duration less than min_duration