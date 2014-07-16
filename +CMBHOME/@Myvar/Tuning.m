function [tuningCurve dim occ] = Tuning(self,variableName,root,varargin)

% Parameters: 
%   dim - Vector of edges for each for the bins
%
%   numBins- If dim is not defined then creates numBins from the minimum of
%            var to max of var (linearly spaced)
%
%   continuizeEpochs - (0) If 1 then take data from all root.epochs. If 0 
%                      then returns output arguments in cell aray, with cell
%                      for each epoch.  

p = inputParser;
p.addParamValue('dim', [], @isnumeric);
p.addParamValue('numBins', 20, @isnumeric);
p.addParamValue('continuizeEpochs', 0, @(c) numel(c)==1 && (c==1 || c==0));

p.parse(varargin{:});
continuizeEpochs = p.Results.continuizeEpochs;


if ~isempty(varargin)
    [occ whichBin dim] = self.Occupancy(variableName,varargin);
else
    [occ whichBin dim] = self.Occupancy(variableName);
end

cel_var = root.cel_myvar2.get(variableName);

if continuizeEpochs
   cel_var = {CMBHOME.Utils.ContinuizeEpochs(cel_var)}; 
end


if (size(occ,2) == 1)
    cel_var = cel_var{1};
    spike_binned = histc(cel_var,dim);
    tuningCurve = spike_binned(1:end-1)./occ;
else
    for i = 1:length(cel_var)
        spike_binned = histc(cel_var{i},dim);
        tuningCurve(:,i) = spike_binned(1:end-1)./occ(:,i);
    end
end
    
end