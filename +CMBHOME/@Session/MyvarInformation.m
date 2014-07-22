function information = MyvarInformation(self, cel, varargin)
% spatial_information = root.SpatialInformation(cel);
%
% Computes the information of a cell about the "myvar" parameter
%
% Will return the information score in bits/spike for
% tetrode cel(1), cell cel(2). Can return continuized epochs, or vectorizes 
% multiple epochs. Note that if dim  
%
% See parameters below.
%
% See Cacucci et al 2007 Methods
%
% OPTIONAL PARAMETERS
%   
%   dim                 vector of bin edges along dimension
%   std_smooth_kernel   STD of the gaussian kernel to smooth the rate map  
%   binside             The length bin when calculating the rate map.
%                       Defaults to 1/20 of the range
% andrew 14 mat 2010

p = inputParser;

p.addRequired('self')
p.addRequired('cel', @isnumeric)
p.addParamValue('dim', [], @isnumeric);
p.addParamValue('std_smooth_kernel', 0, @isnumeric);
p.addParamValue('binside', NaN, @isnumeric)

p.parse(self, cel, varargin{:});

self = p.Results.self;
cel = p.Results.cel;
dim = p.Results.dim;
std_smooth_kernel = p.Results.std_smooth_kernel;
binside = p.Results.binside;

if isnan(binside)
    binside=range(cell2mat2(self.myvar))/20;
end

if isempty(dim)
    dim=min(cell2mat2(self.myvar)):binside:max(cell2mat2(self.myvar));
end

%%
occupancy = histc(cell2mat2(self.myvar),dim)/self.fs_video;
rate_map = histc(cell2mat2(self.spk.myvar),dim)./occupancy;

if std_smooth_kernel>0
    halfWidth = floor(length(rate_map)/2);
    windowWidth = halfWidth*2;
    gaussFilter = gausswin(windowWidth,1/(std_smooth_kernel/binside)*8);
    gaussFilter = gaussFilter / sum(gaussFilter); % Normalize.
    rate_map = conv(rate_map, gaussFilter);
    rate_map=rate_map(hafWidth:end-halfWidth);
end

rate_map=rate_map(occupancy~=0);
occupancy=occupancy(occupancy~=0);

%%

occupancy = occupancy/sum(occupancy);% normalize to probability

F = mean(rate_map);

information =nansum( occupancy .* rate_map .* log2(rate_map ./ repmat(F, size(rate_map, 1), size(rate_map, 2)) )) / F;

end