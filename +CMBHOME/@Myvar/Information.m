function [var_information] = Information(self,variableName,root,varargin)
% var_information = root.myvar2.Information('x',root,...);
%
% Calcuates how much information can be extracted about variable
% 'variableName' from the firing rate of root.cel.  (bits/spike). See
% CMBHOME.Session.SpatialInformation for inspiration.
%
% Inputs:
%   variableName - String of which variable to analyze.
%   root - CMBHOME.Session object. 
%
% Parameters:
%   dim - Vector of edges for each for the bins
%
%   numBins - If dim is not defined then creates numBins from the minimum of
%             var to max of var (linearly spaced)
%
%   continuizeEpochs - (0) If 1 then take data from all root.epochs. If 0 
%                      then returns output arguments in cell aray, with cell
%                      for each epoch.  
%
%   occupation_thresh - Minimum number of samples within a variable bin in
%                       order to count. (0)
%
%   n_thresh - Must have at least this many spikes in the epoch.
%              Otherwise returns NaN (50)

% wchapman 20130729

%% Input parsing
p = inputParser;
p.addParamValue('dim', [], @isnumeric);
p.addParamValue('numBins', 20, @isnumeric);
p.addParamValue('continuizeEpochs', 0, @(c) numel(c)==1 && (c==1 || c==0));
p.addParamValue('occupation_thresh', 0, @isnumeric)
p.addParamValue('n_thresh', 50, @isnumeric)  

p.parse(varargin{:});

continuizeEpochs = p.Results.continuizeEpochs;
occupation_thresh = p.Results.occupation_thresh;
n_thresh = p.Results.n_thresh;

%% Calculation
[tuningCurve,~,occupancy] = self.Tuning(variableName,root);

n_spikes = root.cel_ts;

if iscell(n_spikes)
     n_spikes = cellfun(@length, n_spikes);
    if continuizeEpochs, n_spikes = sum(n_spikes); end
else
    n_spikes = length(n_spikes);
end

for i = 1:size(occupancy, 2)
    occupancy(:, i) = occupancy(:, i) / sum(sum(occupancy(:, i))); % normalize to probability
end

tuningCurve(occupancy <= occupation_thresh) = NaN;
F = nanmean(tuningCurve(:));

var_information = nansum( occupancy .* tuningCurve .* log2(tuningCurve ./ repmat(F, size(tuningCurve, 1), size(tuningCurve, 2)) ), 1) ./ F;

var_information = var_information(:);

var_information(n_spikes<n_thresh) = NaN; % where spiking was too low, get rid of spatial information score

end