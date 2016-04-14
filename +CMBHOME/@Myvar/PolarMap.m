function [tuning_curve theta ang_hd mr u2] = PolarMap(self,variableName,root,varargin)
% [tuning_curve theta ang_hd mr] = PolarMap(self,root,'x')
%
% Tuning of the user defined variable as a function of a polar variable.
%
% Inputs:
%   root - The CMBObject that we lie within (CMBHOME.Session)
%
%   variableName - Which variable to look at (string)
%
% Outputs: 
%   tuning_curve - Vector of mean within each angular bin
%
%   theta - Vector of bin edges (radians, [-pi pi])
%
%   ang_hd - Angle of the mean resultant
%
%   mr - Mean resultant length
%
% Parameters:
%   binsize - (6) Number of bins (radians)
%
%   Continuize - (0) If 1 then averages epochs together. If 1, then return
%                cell array with each cell being an epoch's results.
%
%   ifPlot - (0) Whether or not to plot
%
%   axis - (new) If passed in, then plot on the passed in axes
%
%   polarVar - (root.headdir) The time-varying variable to compare var
%              against


% wchapman 20130724
%% Setup & parseinputs
import CMBHOME.Utils.*
p = inputParser;
p.addParamValue('binsize', 2*pi/60, @(x) isnumeric(x)); % 6 degrees default binsize
p.addParamValue('Continuize',  0, @(x) numel(x)==1);
p.addParamValue('ifPlot', 1, @(x) isnumeric(x));
p.addParamValue('axis', [], @(x) isnumeric(x));

p.parse(varargin{:});
binsize = p.Results.binsize;
Continuize = p.Results.Continuize;
ifPlot = p.Results.ifPlot;
axis = p.Results.axis;

%% Get angles and variable
theta =-pi:binsize:pi;
theta = theta(:);

if Continuize && size(root.epoch,1) > 1
    var = ContinuizeEpochs(self.get(variableName));
elseif size(root.epoch,1)==1
    var = self.get(variableName);
elseif size(root.epoch,1)>1
    var = self.get(variableName);
end

%% get occupancy
[counts wbt] = histc(var,theta);
time = counts / root.fs_video;

%% get spike counts:
[countS wbs] = histc(var(root.cel_i{1}),theta);
%% Create Tuning Curve
tuning_curve = countS ./ time;

%{
% If multiple epochs
if iscell(var)
    tuning_curve = NaN(length(theta)-1,length(var));
    
    for i = 1:length(var)
        [~,wb] = histc(polarVar{i},theta);
        for k = 1:length(tuning_curve)
            tuning_curve(k,i) = nanmean(var{i}(wb==k));
        end
    end
    
% If single epoch    
else
    tuning_curve = NaN(length(theta)-1,1);
    [~,wb] = histc(var,theta);
    for k = 1:length(tuning_curve)
        tuning_curve(k) = nanmean(var(wb==k));
    end
end
%}
%% Stats
if size(tuning_curve,2) == 1
    [ang_hd, mr] = GetOtherStats(tuning_curve, theta); 
else
    for i = 1:size(tuning_curve,2)
        [ang_hd(i), mr(i)] = GetOtherStats(tuning_curve(:,i), theta); 
    end
end

u2 =  CMBHOME.Spike.WatsonsU2(var,var(root.cel_i{1}));

%% ifPlot
if size(tuning_curve,2) > 1, ifPlot=0;end %don't allow user to plot if more than one epoch/cell

if ifPlot
    
   if isempty(axis)
       figure
       axis = gca;
   end

   tc = [tuning_curve(:);tuning_curve(1)];
   theta = [theta(:);theta(1)];
   axes(axis)
   
   polar(theta,tc)
   hold on
   plot([0 cos(ang_hd)*max(tc)],[0 sin(ang_hd)*max(tc)],'r','LineWidth',3)
   
   
   
   title(['U^2=' num2str(u2)])
end

end

%% Get otherStats
function [ang_hd, mr] = GetOtherStats(tuning_curve, theta)

    tf = ~isnan(tuning_curve); % where we arent nan

    %theta=theta*unitsratio('rad','deg');
    
    theta = theta(tf); % remove nans
    
    tuning_curve = tuning_curve(tf);

    xs = tuning_curve.*cos(theta); % average 
    ys = tuning_curve.*sin(theta);

    ang_hd = atan2(mean(ys),mean(xs)); % mean direction
    
    mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(tuning_curve); % mean resultant length
    
end