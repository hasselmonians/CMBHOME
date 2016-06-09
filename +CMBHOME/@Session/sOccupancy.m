function [occupancy, xdim, ydim, pos_locs] = sOccupancy(self, xdim, ydim, mergeepochs, dim)
% Spatially scaled occupancy map
%
% Returns occupancy matrix of time spent in bins of dimension 3x3cm^2 be
% default. occupancy = root.Occupancy(xdim, ydim) where will return a matrix of
% dimensions xdim by ydim(:,2)
%
% Remember that root.spatial_scale indicates cm/pixel in your experiment.
%
% Occupancy is only returned for data within root.epoch. If multiple epochs
% are selected, they are all concatinated and one matrix is returned.
% [occupancy, xdim, ydim] = root.Occupancy returns occupancy matrix and x
% and y dimensions
%
% andrew bogaard 17 may 2010

import CMBHOME.Utils.*

if ~exist('mergeepochs', 'var'), mergeepochs = 1; end

if ~exist('dim', 'var'), dim = 3; end % cm width of each bin ( dim^2 cm^2 bins)

if mergeepochs
    self = MergeEpochs(self); 
    x=0;y=0;
    while(range(x)==0 || range(y)==0)
        figure(100);
        plot(self.sx,self.sy);
        close(100);
    x = cell2mat2(self.sx);
    y = cell2mat2(self.sy);
    end
else
    x=0;y=0;
    while(range(x)==0 || range(y)==0)
    x = self.sx;
    y = self.sy;
    end
end

if ( ~exist('xdim', 'var') || ~exist('ydim', 'var') ) || (isempty(xdim) || isempty(ydim))
    if iscell(x)
        xdim = min(vertcat(x{:})):dim:max(vertcat(x{:})); %edges of x and y dimensions

        ydim = min(vertcat(y{:})):dim:max(vertcat(y{:}));
    else
        xdim = min(x):dim:max(x); %edges of x and y dimensions

        ydim = min(y):dim:max(y);
    end
    
end

if ~iscell(x)

    %     occupancy = hist3([x, y], 'Edges', {xdim, ydim}) / self.fs_video;
    [occupancy,~,~,pos_locs] = histcn([x y],[xdim max(xdim+dim)],[ydim max(ydim+dim)]);
    occupancy = occupancy/self.fs_video;
%     keyboard
    
    %%
%     for i=1:size(pos_locs,1)
%        temp(pos_locs(i,2),pos_locs(i,2))=1; 
%     end
%     
%     %%
%     keyboard

    temp = pos_locs;
    
    %
    pos_locs = zeros(size(pos_locs,1),1);
    pos_locs(temp(:,1)>0&temp(:,1)<size(occupancy,1)&temp(:,2)>0&temp(:,2)<size(occupancy,2))=...
        sub2ind(size(occupancy),...
        temp(temp(:,1)>0&temp(:,1)<size(occupancy,1)&temp(:,2)>0&temp(:,2)<size(occupancy,2),1),...
        temp(temp(:,1)>0&temp(:,1)<size(occupancy,1)&temp(:,2)>0&temp(:,2)<size(occupancy,2),2));
    
else

    occupancy = zeros(length(xdim), length(ydim), length(x));
    pos_locs = cell(length(x),1);
    
    for i = 1:length(x)
        
        %         occupancy(:, :, i) = hist3([x{i}, y{i}], 'Edges', {xdim, ydim}) / self.fs_video;
        [occupancy(:,:,i),~,~,pos_locs{i}] = histcn([x{i} y{i}],[xdim max(xdim+dim)],[ydim max(ydim+dim)]);
        occupancy(:,:,i) = occupancy(:,:,i)/self.fs_video;
        pos_locs{i} = sub2ind(size(occupancy),pos_locs{i}(:,1),pos_locs{i}(:,2));
        
    end
    
end