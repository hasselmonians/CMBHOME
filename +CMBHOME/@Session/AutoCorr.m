function [cor, lag, smooth, theta_skip_index] = AutoCorr(self, cel, varargin)
% Returns and optionally plots spiketime autocorrelation
%
% Capable of velocity thresholding, and multiple epoch averaging. If root
% has multiple epochs and average_epochs is set to 1, then the returned
% correlation is a vector. If it is set to 0, then a column array is
% returned.
%
% If velocity thresholding is set, and average_epochs==1, then only epochs
% with running speed that meet speed_thresh are calculated. If multiple
% epochs are set, only running epochs within those epochs are solved for.
% 
% If multiple epochs are set, and velocity thresholding is set, and
% average_epochs==0, then there is a good chance that the returned cor
% array will have as many columns as running epochs within root.epoch.
%
% PARAMS:
%   speed_thresh                [min_speed max_speed] (default [-1 -1] for no speed thresh
%   theta_skip                  calculate and print to plot (if plotting, see below) theta skipping score (0) 
%   t_bin                       (seconds) time binwidth (.020)
%   max_lag                     (seconds) maximum lag (.6) 
%   plot_smoothed_signal        plots black line above bars indicating smoothed signal (1)
%   average_epochs              0 or 1, (1). if 1, returns one single autocorrelation, if 0, 
%                               returns autocorr for each epoch.
%   supress_plot                if 1, does not plot (0)
%   figure_handle               if defined, plots to figure with given
%                               handle
%   std_smooth_kernel           (seconds) std of kernel to smooth xcorr
%                               signal (default is 2*t_bin)
%   unbiased                    0 or 1 (1). if 0, returns biased xcorr
%   center_peak                 0 or 1 (0). if 1, keeps center peak, if 0,
%                               sets center peak to 0
%   
% ALGORITHM: tbc
%
% andrew oct 20 2010
% v 1.1 jan 12 2011.    removed user_kalman_vel param, because now the object
%                       defaults to the user defined root.b_vel, if it exists
%
% [cor, lag, smooth, theta_skip_index] = root.AutoCorr(cel, varargin)

p = inputParser;

p.addRequired('self')
p.addRequired('cel', @isnumeric)
p.addParamValue('speed_thresh', [-1 -1], @(c) (numel(c)==2 && diff(c)>=0));
p.addParamValue('theta_skip', 0, @(c) (c==1 || c==0));
p.addParamValue('t_bin', .020, @(c) numel(c)==1);
p.addParamValue('max_lag', .6, @(c) numel(c)==1);
p.addParamValue('plot_smoothed_signal', 1, @(c) (c==1 || c==0));
p.addParamValue('supress_plot', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('figure_handle', '', @(c) numel(c)==1);
p.addParamValue('std_smooth_kernel', [], @isnumeric);
p.addParamValue('average_epochs', 1, @(c) (c==1 || c==0));
p.addParamValue('unbiased', 1, @(c) (c==1 || c==0));
p.addParamValue('center_peak', 0, @(c) (c==1 || c==0));

p.parse(self, cel, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
speed_thresh = p.Results.speed_thresh;
theta_skip = p.Results.theta_skip;
t_bin = p.Results.t_bin;
max_lag = p.Results.max_lag;
plot_smoothed_signal = p.Results.plot_smoothed_signal;
supress_plot = p.Results.supress_plot;
figure_handle = p.Results.figure_handle;
std_smooth_kernel = p.Results.std_smooth_kernel;
average_epochs = p.Results.average_epochs;
unbiased = p.Results.unbiased;
center_peak = p.Results.center_peak;

theta_skip_index = [];

if isempty(std_smooth_kernel), std_smooth_kernel = t_bin*2; end

import CMBHOME.Utils.*

% set epochs to appropriate running epochs

[self, success] = RunningEpochs(self, speed_thresh(1), speed_thresh(2), max_lag);

if ~success
    text(.2, .4, 'No epochs meet this running speed requirement'); 
    axis off; 
    return; 
end

% Calculate xcorrelations

[cor, lag] = self.spk_xcorr(cel, max_lag, t_bin, average_epochs); % average autocorrelation where running speed is above min_vel

%if ~unbiased, cor = Biased(cor); end

[cor, lag, smooth] = SmoothCor(cor, lag, t_bin, std_smooth_kernel); % solve for smoothed signal, and chop off center peak and negative half of autocorrelation

if theta_skip, theta_skip_index = self.ThetaSkipping(cel, 1); end % calc theta skip, if requested

if ~supress_plot, PlotIt(cor, lag, smooth, plot_smoothed_signal, theta_skip_index, figure_handle); end
    
end    

function PlotIt(cor, lag, smooth, plot_smoothed_signal, theta_skip_index, figure_handle)

    if ~isempty(figure_handle), figure(figure_handle); end
    
    if isempty(cor)
        text(.2, .4, 'No auto correlation for these epochs.'); 
        axis off; 
        return;
    end
    
    if isvector(cor)
    
        bar(lag, cor,1, 'FaceColor', [.2 .2 1], 'EdgeColor', [ 0 0 .75]), hold on;
        
        if plot_smoothed_signal, line(lag, smooth, 'Color', 'k', 'linewidth', 1.5); end
        
        if lag(end)>0
            xlim([0 lag(end)]);
        end

        if max(smooth)>0
            ylim([0 max(smooth)*1.1]);
        end
        
    else % multiple autocorrelations to plot
    
        import CMBHOME.Utils.*
                
        line(lag, cor), hold on;

        if lag(end, end) > 0
            xlim([0 lag(end, end)]);
        end
        
    end
    
    ys = ylim;
    xs = xlim;

    if ~isempty(theta_skip_index), text(xs(2)*.6, ys(2)*.9, ['Theta Skipping: ' num2str(theta_skip_index)], 'FontSize', 8); end

    %title('Auto Correlation of Spike Train');
    
    %xlabel('Lag (Seconds)');
    
end

function cor = Biased(cor)
% re-bias the unbiased xcorr

    lags = (-floor(length(cor)/2):floor(length(cor)/2));

    scale = (ceil(length(cor)/2)-abs(lags));
	scale(scale<=0)=1; % avoid divide by zero, when correlation is zero

    cor = cor(:).*scale(:);

end

function [cor, lag, smooth] = SmoothCor(cor, lag, t_bin, std_smooth_kernel)
    
    import CMBHOME.Utils.*

    std_smooth_kernel = .005;
    
    kernel = pdf('normal', -std_smooth_kernel*10/t_bin:std_smooth_kernel*10/t_bin, 0, std_smooth_kernel/t_bin); % convolve with gaussian

    if isempty(cor), smooth = []; return; % if no signals, then send it back
    
    else
        
        smooth = zeros(size(cor));
        
        for i = 1:size(cor,2)  

            smooth(:,i) = ndnanfilter(cor(:,i), kernel(:), []);
        
        end
        
        lag = lag(round(end/2):end,:);
        
        cor = cor(round(end/2):end,:);
        
        smooth = smooth(round(end/2):end, :);

    end 
    
end

function [self, success] = RunningEpochs(self, lower_threshold, upper_threshold, speed_dur)

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