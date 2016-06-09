function [self, success] = RunningEpochs(self, lower_threshold, upper_threshold, speed_dur)
    % Sets lower and upper thresholds for running speed & sets epochs accordingy
    import CMBHOME.Utils.*
    
    success = 0;

    if upper_threshold <= 0, success = 1; return; end
    
    if iscell(self.vel) % scale velocity to cm/sec
        
        vel = cellfun(@(c) c*self.spatial_scale, self.vel, 'Unif', 0);
        
    else
        
        vel = self.vel * self.spatial_scale;
        
    end
        
    ind = ThresholdBandDetect(vel, lower_threshold, upper_threshold, 1, speed_dur*self.fs_video);   % returns all indeces for epochs longer than .5 seconds that meet threshold
                                                                                                    % from +CMBHOME/+Utils
    if isempty(ind)
        disp('No epochs that meet requirements for running speed');
        return;
    end
    
    if iscell(ind)
        
        epoch = [0 0];
        
        for i = 1:length(ind), epoch = cat(1, epoch, self.ts{i}(ind{i})); end
        
        epoch(1,:) = [];
        
    else

        epoch = self.ts(ind); % set epochs to valid running epochs
        
    end
    
    if size(epoch, 2)~=2, epoch = epoch'; end

    success = 1;
    
    self.epoch = epoch;
    
end   