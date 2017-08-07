function [self, success] = regionEpochs(self, region, region_dur)
    % Sets lower and upper thresholds for running speed & sets epochs accordingy
    import CMBHOME.Utils.*
    
    success = 0;
    
    x = self.b_x;
    y = self.b_y;
    
    inds_x = x>=region(1) & x<=region(2);
    inds_y = y>=region(3) & y<=region(4);
    inds = inds_x & inds_y;
        
    ind = ThresholdBandDetect(inds, 0, 1, 1, region_dur*self.fs_video);   % returns all indeces for epochs longer than .5 seconds that meet threshold
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
    
    epoch = FilterEpochs(epoch, region_dur, 5/self.fs_video);
    
    self.epoch = epoch;
    
end   