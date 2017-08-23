function [self, success] = regionEpochs(self, region, region_dur)
    % Sets lower and upper thresholds for running speed & sets epochs accordingy
    import CMBHOME.Utils.*
    
    success = 0;
    
    x = self.b_x;
    y = self.b_y;
    
    inds_x = x>=region(1) & x<=region(2);
    inds_y = y>=region(3) & y<=region(4);
    inds = inds_x & inds_y;
   
    [ind, tf] = CMBHOME.Utils.OverThresholdDetect(inds, 0.5, 1, 1);

    %{
    indd = zeros(size(x));
    for i = 1:size(ind,1)
        indd(ind(i,1):ind(i,2)) = 1;
    end
    indd = indd==1;
    %figure
    clf
    subplot(2,2,1)
    plot(x,y,'.')
    xlim([nanmin(x) nanmax(x)]); ylim([nanmin(y) nanmax(y)])
    subplot(2,2,2)
    plot(x(inds_x), y(inds_x),'.');
    xlim([nanmin(x) nanmax(x)]); ylim([nanmin(y) nanmax(y)])
    subplot(2,2,3)
    plot(x(inds_y), y(inds_y),'.');
    xlim([nanmin(x) nanmax(x)]); ylim([nanmin(y) nanmax(y)])
    subplot(2,2,4)
    plot(x(inds),y(inds),'.');
    xlim([nanmin(x) nanmax(x)]); ylim([nanmin(y) nanmax(y)])
    hold on
    plot(x(indd),y(indd),'r.')
    xlim([nanmin(x) nanmax(x)]); ylim([nanmin(y) nanmax(y)])
    %}
    %%
    
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