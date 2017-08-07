function self = FixPos(self, max_allowed_flips)
% Fixes large jumps in position and interpolates missing values
%
% Takes all 0,0s and large jumps in position (greater than
% jitter_threshold=15 pixels/sample) that persist for less than
% max_allowed_flips (default = 5) and linearly interpolates the missing
% data. Smooths conservatively afterward, as well (convolution with a gaussian, standard
% deviation = 2 samples).
%
% root = root.FixPos;

import CMBHOME.Utils.*

if self.raw_pos~=1
    disp('It appears the tracking data has already been affected. Reset root.raw_pos=1 to override.');
    return;
end
    
warning('off', 'MATLAB:interp1:NaNinY');

    jitter_threshold = 10/self.spatial_scale; % pixels in 10 cm; one-sample change in distance that qualifies as bad vector
    
    x = self.b_x;
    y = self.b_y;
    
    ts = self.b_ts;
    
    bads = (x==0 | y==0);
    
    x(bads) = NaN;
    y(bads) = NaN;
      
    if ~exist('max_allowed_flips', 'var')
        max_allowed_flips = 5; % samples
    end
    
    flips = findOnsetsAndOffsets(isnan(x));
    
    flips(:,2) = flips(:,2)+1;
    
    flips = cat(1, 1, flips(:), find([0; sqrt(diff(x).^2 + diff(y).^2)]>jitter_threshold), length(x));
    
    flips = sort(unique(flips)); % indeces of NaNs or jumps
    
    flips = [flips(1:end-1), flips(2:end)];  % epochs formation
    
    flips(:,2) = flips(:,2)-1; % adjust for diff shift
    
    flips(flips(:,2)-flips(:,1)>max_allowed_flips-1,:) = [];
    
    flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps
   
    flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices
    
    x([flips{:}]) = []; % remove samples in ts and x
    y([flips{:}]) = []; % remove samples in ts and x
    ts([flips{:}]) = [];
    
    x = interp1(ts, x, self.b_ts);
    y = interp1(ts, y, self.b_ts);

    x = ndnanfilter(x, normpdf(-6:6, 0, 2)', [], 1, {}, {}, 1); % conv with gaussian and ignore NaNs
    y = ndnanfilter(y, normpdf(-3:3, 0, 2)', [], 1, {}, {}, 1);
    
    self.raw_pos = 0 ;
    self.b_x = x;
    self.b_y = y;

    self.b_x = self.b_x - min(self.b_x);
    self.b_y = self.b_y - min(self.b_y);
    
end
