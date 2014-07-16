function [varSpeed speedDim] = speedModulation(self,variableName,root,varargin)
% [varSpeed speedDim] = root.myvar2.speedModulation(root,'x','ifplot',1)
%
% Calculates the mean of the variable within a certain speed range
%
% Outputs:
%   varSpeed -- The mean within each epoch (vector / cell)
%   speedDim -- The edges of the speed bins (vector)
%
% Inputs:
%   root: The CMBHOME.Session object
%   variableName: Which variable within root.myvar2 to look at (string)
%
% Parameters:
%   speedDim:
%   numBins:
%   continuize_epochs: 
%   ifplot: (0) Whether or not to plot the results
%   plotAxis: (new axis) If passed and ifplot=1 then plots on the passed
%             axes

% wchapman 20130729

%% Parse inputs
p = inputParser;
p.addParamValue('speedDim',[])
p.addParamValue('numBins',[])
p.addParamValue('continuize_epochs', 0, @(c) numel(c)==1 && (c==1 || c==0))
p.addParamValue('ifplot', 0, @(c) numel(c)==1 && (c==1 || c==0))
p.addParamValue('plotAxis',[],@isnumeric)

p.parse(varargin{:});

speedDim = p.Results.speedDim;
numBins = p.Results.numBins;
continuize_epochs = p.Results.continuize_epochs;
ifplot = p.Results.ifplot;
plotAxis = p.Results.plotAxis;

%% Get the variable
var = self.get(variableName);
speed = root.vel;

if iscell(var)
    if continuize_epochs
        var = CMBHOME.Utils.ContinuizeEpochs(var);
        speed = CMBHOME.Utils.ContinuizeEpochs(speed);
    end
end

isCirc = self.isCirc_get(variableName);

%% Get the speed bins
if iscell(speed), speed2 = CMBHOME.Utils.ContinuizeEpochs(speed); else speed2 = speed;end

if isempty(speedDim)
    if ~isempty(numBins)
        speedDim = linspace(0,max(speed2),numBins);
    else
        speedDim = 0:root.spatial_scale^-1:max(speed2);
        numBins = length(speedDim)-1;
    end
end

%% Analysis

if ~iscell(var)
    [~,whichBin] = histc(speed,speedDim);
    varSpeed = NaN(max(whichBin)-1,1);
    for i = 1:numBins
        inds = whichBin == i;
        varSpeed(i) = nanmean(var(inds));
    end
else
    varSpeed = cell(length(var),1);
    for k = 1:length(var)
        [~,whichBin] = histc(speed{k},speedDim);
        varSpeed{k} = NaN(max(whichBin)-1,1);
        for i = 1:numBins
            inds = whichBin == i;
            varSpeed{k}(i) = nanmean(var{k}(inds));
        end
    end
end

speedDim = speedDim(1:end-1);
speedDim = speedDim(:);

%% Ifplot
if ifplot && ~iscell(var)
    if isempty(plotAxis)
        figure
        plotAxis = gca;
    end
    axes(plotAxis)
    plot(speedDim,varSpeed,'.')
    
    xlabel('Running Speed (pixels/second)')
    ylabel(variableName)
    
    if iscell(root.name)
        title(root.name{end})
    else
        title(root.name);
    end
    
end

end