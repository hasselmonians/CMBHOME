function [ spkpos,spk_i ] = spk_pos( vid_ts,pos,spk_ts )
%SPK_POS Finds spike positions
%   Finds the spike positions given tracking data in a N-dimensional state
%   space
%
%   SPKPOS = SPK_POS(VID_TS,POS,SPK_TS);
%
%   ARGUMENTS
%   * VID_TS: A vector of the timestamps for the state data
%   * POS: The position at every timestamp. A NXM matrix, where N is the
%   number of time samples in VID_TS and M is the number of dimensions of
%   the space
%   * SPK_TS: The times of the spikes
%
%   RETURNS
%   * SPKPOS: The state of the animal at the time of each spike
%   * SPK_I: A vector such that SPKPOS = POS(SPK_I,:)
%
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com)
i = 1:numel(vid_ts);
spk_i = interp1(vid_ts,i,spk_ts,'nearest','extrap');
spkpos = pos(spk_i,:);

end

