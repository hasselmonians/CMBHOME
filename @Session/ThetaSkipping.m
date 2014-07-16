function R = ThetaSkipping(self, cel, supress_plot)
%
% R = root.ThetaSkipping(cel);
%
% Calculates the ratio of the difference between 1st and 2nd peaks in
% autocorrelogram to the max of both peaks is R
%
% See PAGE 996 in Sachin S. Deshmukh, D. Yoganarasimha, Horatiu Voicu and James J. Knierim
% J Neurophysiol 104:994-1006, 2010
%
% andrew abogaard 10 sept 2010

    import CMBHOME.Utils.*
    
    if ~exist('supress_plot', 'var'), supress_plot = 0; end
            
    % autocorrelation parameters
    max_lag_t = .6; % seconds
    t_bin = .005; % seconds
    speed_thresh = [5 100000];
    speed_dur = max_lag_t;

    % calculate autocorrelation
    import CMBHOME.Utils.*

    %kernel = pdf('normal', -100:100, 0, .05/t_bin);
    
    % set epochs to appropriate running epochs

    [cor, lag, epochs] = self.spk_xcorr(cel, max_lag_t, t_bin, 1);

    cor(lag==0) = [];

    lags = lag(round(end/2):end);
    %acs = conv(cor(round(end/2):end), kernel, 'same');

        % filter signal between 1 and 10 hz
        Wn_theta = [1/(t_bin^-1/2) 10/(t_bin^-1/2)]; % normalized by the nyquist frequency

        [btheta,atheta] = butter(3,Wn_theta);
        
        acs = filtfilt(btheta,atheta,cor);
        
        acs = acs(round(end/2):end)+mean(cor);
    
    [amp, ind] = extrema(acs);
    
    l = lags(ind);
    
    peak1 = max(amp(l>=.1 & l<=.2)); % first peak
    
    peak2 = amp(find(l>.2, 1, 'first')); % second peak
    
    if ~isempty(peak1) && ~isempty(peak2)
    	R = (peak2-peak1) / max([peak1 peak2]); % theta skipping index
    else
        R = []; 
        disp('No peaks in theta range of autocorrelogram');
    end
    
end