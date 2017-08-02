function [diffData2] = circDiff(data,dim,type)
% function [diffData] = circDiff(data,dim,type)
%
% computes the difference between each point accounting for the
% fact that the values are circular.
%  INPUTS:
%  data - circular data type vector or matrix
%  dim  - dimension to operate on, default=1
%  type - 'deg' or 'rad' to specify if 2pi or 360 is max value

data = data - min(data(:));
import CMBHOME.Axona.*
if ~exist('type','var')
	if abs(360-max(data(:))) > abs(2*pi-max(data(:)))
		% data are probably in radians
		cycle = 2*pi;
	else
		% data are probably in degrees
		cycle = 360;
	end
else
	switch(type)
		case 'deg'
			cycle = 360;
		case 'rad'
			cycle = 2*pi;
		otherwise
			error('Only use deg or rad for type');
	end
end
	
if ~exist('dim','var'), dim = 1; end;

diffData = diff(data,[],dim);
dataDims = ndims(diffData);
diffData = cat(dataDims+1,diffData, diffData+cycle, diffData-cycle); 

diffData2 = minMag(diffData,dataDims+1);


