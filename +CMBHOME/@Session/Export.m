function [dat, beh, lfp, spk] = Export(self, fname)
    % Saves all important (raw) data as a simple struct. Useful for exporting
    % to other frameworks / languages.

    %%
    beh_x = self.b_x;
    beh_y = self.b_y;
    beh_hd = self.b_headdir;
    beh_spatial_scale = self.spatial_scale;
    beh_vel = self.b_vel;
    beh_ts = self.b_ts;

    %%
    for i = 1:length(self.b_lfp)
        if ~isempty(self.b_lfp(i).signal)
            lfp_signal(i,:) = self.b_lfp(i).signal;
            lfp_ts(i,:) = self.b_lfp(i).ts;
            lfp_fs(i) = self.b_lfp(i).fs;
            lfp_chan(i) = i;
        end
    end


    %%
    for i = 1:size(self.spike,1)
        for k = 1:size(self.spike,2)
            if ~isempty(self.spike(i,k).i)
                spk_ts(i,k) = {self.spike(i,k).ts};
            end     
        end
    end
    
    %%
    if ~exist('fname','var')
        [pn,pd] = uiputfile('*.mat');
        fname = [pd, pn];
    end
    
    %%
    if length(isempty(self.b_lfp)) && isempty(self.b_lfp(1).signal)
        save(fname,'spk_ts',...
            'beh_hd','beh_spatial_scale','beh_ts','beh_vel','beh_x','beh_y',...
            'spk_ts')
    else
        save(fname,'spk_ts',...
            'beh_hd','beh_spatial_scale','beh_ts','beh_vel','beh_x','beh_y',...
            'spk_ts', ...
            'lfp_chan','lfp_ts','lfp_signal','lfp_ts','lfp_fs')
    end
    
end