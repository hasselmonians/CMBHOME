function [occ whichBin dim] = Occupancy(self,variableName,varargin)
% [occ whichBin dim] = root.myvar2.Occupancy('x');
%
% Returns occupancy of root.myvar in each of the bins. 
%
% Inputs: 
%   variableName- string of which variable to bin
%
% Outputs:
%   occ - Number of samples in each bin
%
%   whichBin - Tells which bin each of the samples of var is in
%
%   dim - Edges for each of the bins
%
% Parameters: 
%   dim - Vector of edges for each for the bins
%
%   numBins- If dim is not defined then creates numBins from the minimum of
%            var to max of var (linearly spaced)
%
%   continuizeEpochs - (0) If 1 then take data from all root.epochs. If 0 
%                      then returns output arguments in cell aray, with cell
%                      for each epoch.  

% wchapman 20130729

p = inputParser;
p.addParamValue('dim', [], @isnumeric);
p.addParamValue('numBins', 20, @isnumeric);
p.addParamValue('continuizeEpochs', 0, @(c) numel(c)==1 && (c==1 || c==0));

p.parse(varargin{:});
dim = p.Results.dim;
numBins = p.Results.numBins;
continuizeEpochs = p.Results.continuizeEpochs;

if(self.isCirc_get(variableName))
   dim = -pi:(2*pi/numBins):pi;
end


var = self.get(variableName);

if continuizeEpochs
    var = CMBHOME.Utils.ContinuizeEpochs(var);
end

if ~iscell(var)
    if isempty(dim)
        dim = linspace(min(var),max(var),numBins+1);
    end
    [occ whichBin] = histc(var,dim);
    occ = occ(1:end-1);
else
    var2 = CMBHOME.Utils.ContinuizeEpochs(var);
    if isempty(dim)
        dim = linspace(min(var2),max(var2),numBins+1);
    end
    occ = NaN(length(dim)-1,length(var));
    whichBin = cell(length(var),1);
    for i = 1:length(var)
        [occt whichBin{i}] = histc(var{i},dim);
        occ(:,i) = occt(1:end-1);
    end
end

dim = dim(:);

end