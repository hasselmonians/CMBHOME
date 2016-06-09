function [f_intrins, theta_index, modelstruct] = IntrinsicFrequency2(self, cel, varargin)
% Intrinsic frequency of a unit
%
% ARGUMENTS
%
%   cel                     [tetrode index, cell index] in Session object root
%   params                  (see below)
%
% RETURNS
%
%   f_intrins               scalar, frequency of fit line
%   theta_index             scalar, measure of theta-ness
%   modelstruct             struct, with fields:
%
%                             model
%                             cor
%                             lag
%                             theta_index
%                             f_intrins
%                             theta_skipping
%
% PARAMS:
%
%   speed_thresh                [min_speed max_speed] (default [-1 -1] for no speed thresh
%   theta_skip                  calculate and print to plot (if plotting, see below) theta skipping score (0) 
%   t_bin                       (seconds) time binwidth (.020)
%   max_lag                     (seconds) maximum lag (.6) 
%   average_epochs              0 or 1, (1). if 1, returns one single autocorrelation, if 0, 
%                               returns autocorr for each epoch. (averages xcorr signal 
%                               weighted by duration of epoch)
%   supress_plot                if 1, does not plot (0)
%   figure_handle               if defined, plots to figure with given
%                               handle
%
% THE ALGORITHM
%
% See Royer et al 2010 J. Neuroscience, Distinct representations and theta dynamics in dorsal and ventral hippocampus. pg 1781 and Figure 6G
%
% Calculate unbiased xcorr of spike trains that 
% occur during epochs of running longer than .5 seconds. Fit the xcorr
% signal to the equation in Royer et al 2010, return the omega value of the
% fit, as well as the ratio of a/b as thetaness score
%
% andrew bogaard 25 oct 2010
% v 1.1 jan 12 2011.    removed user_kalman_vel param, because now the object
%                       defaults to the user defined root.b_vel, if it
%                       exists
% v 1.2 july 4 2011     updated to measure theta skipping!
%
% (1) [f_intrins, theta_index] = root.IntrinsicFrequency2(cel);
% (2) [f_intrins, theta_index, modelstruct] = root.IntrinsicFrequency2(cel, params)

p = inputParser;

p.addRequired('self')
p.addRequired('cel', @isnumeric)
p.addParamValue('speed_thresh', [-1 -1], @(c) (numel(c)==2 && diff(c)>=0));
p.addParamValue('theta_skip', 0, @(c) (c==1 || c==0));
p.addParamValue('t_bin', .01, @(c) numel(c)==1);
p.addParamValue('max_lag', .7, @(c) numel(c)==1);
p.addParamValue('supress_plot', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('figure_handle', '', @(c) numel(c)==1);
p.addParamValue('average_epochs', 1, @(c) (c==1 || c==0));
p.addParamValue('theta_skipping', 0, @(c) (c==1 || c==0));

p.parse(self, cel, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
speed_thresh = p.Results.speed_thresh;
theta_skip = p.Results.theta_skip;
t_bin = p.Results.t_bin;
max_lag = p.Results.max_lag;
supress_plot = p.Results.supress_plot;
figure_handle = p.Results.figure_handle;
average_epochs = p.Results.average_epochs;
theta_skipping = p.Results.theta_skipping;

%[cor, lag] = self.AutoCorr(cel, 'speed_thresh', speed_thresh, 't_bin', t_bin, 'max_lag', max_lag,'supress_plot',1);

if t_bin / mod(max_lag, t_bin) ~= 2 % set lags so it is 'even' (odd number of coefficients and zero centered')

    max_lag = t_bin*floor(max_lag/t_bin)+.5*t_bin;
    
end
    
[cor, lag] = self.spk_xcorr(cel, max_lag, t_bin, 1, 'prob');

if theta_skipping % new model
    
    f = fittype('[a*(cos(w*x)+1)+a2*(cos(w*x/2)+1)+b]*exp(-abs(x)/tau)+c*exp(-x^2/tau2^2)',...
            'independent', 'x', ...
            'coefficients', {'a' 'a2' 'b' 'c' 'tau' 'tau2' 'w'}); % make model
    
    options = fitoptions(f); 
    
    options.Lower = [0 0 0 -inf 0 0  5*pi*2];
    options.Upper = [1 1 1 inf 5 .050 9*pi*2];
    options.MaxIter = 10^3;
    options.MaxFunEvals = 10^4;
    
    options.StartPoint = [  max(cor(lag>.1 & lag<.150))-min(cor(lag>.1 & lag<.150)),... % a
                            0,... % a2
                            mean(cor(lag>-.1&lag<.1)),... %b
                            cor(round(end/2)),... %c
                            .1,... % tau
                            .015,... % tau2
                            2*pi*8]; % w

    [model, gof, info] = fit(lag(lag~=0), cor(lag~=0), f, options);

    peak_ts = 2*model.a+2*model.a2+model.b; % peak at theta skipping peak (3rd peak)
    
    peak_t = 2*model.a+model.b; % peak at theta skipping trough (2nd peak)
    
    %SL(i).thetaskip2 = model.a2 / model.a; looked pretty good
    
    theta_skipping = (peak_ts-peak_t) / max(peak_t, peak_ts);
    
    theta_index = (model.a+model.a2) / mean(cor);
    
else % original analysis
    
    cor = cor / max(cor(lag>.1 & lag<.150)); % normalize to peak between 100 and 150ms

    cor(cor>1) = 1; % clip anything greater than 1
    
    f = fittype('[a*(cos(w*x)+1)+b]*exp(-abs(x)/tau)+c*exp(-x^2/tau2^2)'); % make model
    
    options = fitoptions(f); 

    options.Lower = [0 0 0 0 0  5*pi*2];
    options.Upper = [1 1 1 5 .050 9*pi*2];
    options.MaxIter = 10^3;
    options.MaxFunEvals = 10^4;
    %                      
    options.StartPoint = [  max(cor(lag>.1 & lag<.150))-min(cor(lag>.1 & lag<.150)),... % a
                            mean(cor(lag>-.1&lag<.1)),... %b
                            cor(round(end/2)),... %c
                            .1,... % tau
                            .015,... % tau2
                            2*pi*8];  % w
    
    [model, gof, info] = fit(lag(lag~=0), cor(lag~=0), f, options);
    
    theta_skipping = NaN; % didnt measure for theta skipping

    theta_index = model.a / mean(cor);
    
end

if info.exitflag<=0 % no convergence
   
    disp('No convergence in IntrinsicFrequency2');
    
    theta_index = NaN;
    f_intrins = NaN;
    modelstruct.model = model;
    modelstruct.gof = gof;
    modelstruct.info = info;
    modelstruct.cor = cor;
    modelstruct.lag = lag;
    modelstruct.theta_index = NaN;
    modelstruct.f_intrins = NaN;
    modelstruct.theta_skipping = NaN;
    
else

    f_intrins = model.w/(2*pi);

    modelstruct.model = model;
    modelstruct.gof = gof;
    modelstruct.info = info;
    modelstruct.cor = cor;
    modelstruct.lag = lag;
    modelstruct.theta_index = theta_index;
    modelstruct.f_intrins = f_intrins;
    modelstruct.theta_skipping = theta_skipping;

end

if ~supress_plot
    
    PlotIt(model, cor, lag, figure_handle);
    
end
end

function PlotIt(model, cor, lag, figure_handle)
% plots root.IntrinsicFrequency2 results

if ~isempty(figure_handle), figure(figure_handle); end

bar(lag, cor, 'edgecolor', 'k', 'facecolor', 'k', 'barwidth', 1), 

xlim([min(lag) max(lag)]), hold on

ca = plot(model, 'r-');
set(ca, 'LineWidth', 2.5)
end
