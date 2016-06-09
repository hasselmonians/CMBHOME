function x = RemoveJumpsAndSmooth(x, max_allowed_flips, jitter_threshold, smooth_std,filter_def)
%  Removes large jumps in data, smoothes the removed regions
%   ARGUMENTS
%       x                   a vector of something. NaNs may be used in place of
%                           known bad samples
%       max_allowed_flips   (samples) How many consecutive samples can be
%                           corrected? default = 5
%       jitter_threshold    (units of x) What intersample interval
%                           constitutes a jump?
%                           default = 15
%       smooth_std          (samples) STD of gaussian smoothing kernel.
%                           default = 2
%
% This function takes a vector, x, and looks for jumps greather than
% jitter_threshold that last for less than max_allowed_flips in samples. We
% then linearly interpolate the missed samples and smooth by convolution
% with a gaussion of standard deviation smooth_std bins.
%
% x = RemoveJumpsAndSmooth(x, max_flip_length, jitter_threshold, smooth_std);

    import CMBHOME.Utils.*

    x = x(:);
    
    ts = 0:length(x)-1;
    
    ts2 = ts;
    
    if ~exist('jitter_threshold', 'var'), jitter_threshold = 15; end % pixels; one-sample change in distance that qualifies as bad vector
    
    if ~exist('max_allowed_flips', 'var'), max_allowed_flips = 5; end % samples
    
    if ~exist('smooth_std', 'var'), smooth_std = 2; end % samples

    if ~exist('filter_def','var'), filter_def = -6:6;end
    
    % there are two cases for which are the flips, and which is the real data. we will see how many
    % samples fall into each case, and then pick the lesser. then we can
    % cut out any jumps that are too long
    
    flips = findOnsetsAndOffsets(isnan(x));
    
    flips(:,2) = flips(:,2)+1;
    
    flips = [flips(:); find(abs(diff(x))>jitter_threshold)+1]; % all flip points

    flips = sort(unique(flips)); % indeces of NaNs or jumps
    
    flips = [flips(1:end-1), flips(2:end)];  % epochs formation
    
    flips(:,2) = flips(:,2)-1; % adjust for diff shift
    
    flips(flips(:,2)-flips(:,1)>max_allowed_flips-1,:) = []; % remove those greater than max_allowed_flips
    
    flips = mat2cell(flips, ones(size(flips,1),1),2); % convert to pairs corresponding to steps
   
    flips = cellfun(@(c) c(1):c(2), flips, 'unif', 0); % convert to indices
    
    x([flips{:}]) = []; % remove samples in ts and x
    ts([flips{:}]) = [];
    
    x = interp1(ts, x, ts2);
    
    
    % gaussian_filter = pdf('Normal',-5*smooth_std:5*smooth_std,0,smooth_std);
    if smooth_std>0
        x = ndnanfilter(x, normpdf(filter_def, 0, smooth_std)', [], 1, {}, {}, 1); % convolve with gaussian and ignore NaNs
    end
    % x = conv(x,gaussian_filter, 'same'); % smooth by convolution with a gaussian
    
end
