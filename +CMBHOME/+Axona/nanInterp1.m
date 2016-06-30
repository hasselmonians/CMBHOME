function [ varargout ] = nanInterp1( varargin )
%function [ varargout ] = nanInterp1( varargin )

%   Detailed explanation goes here

import CMBHOME.Utils.*

varargin(cellfun(@isnumeric,varargin)) = cellfun(@nanInterp,varargin(cellfun(@isnumeric,varargin)),'UniformOutput',false);
varargout = cell(1,nargout);

[~,temp]=unique(varargin{1});
varargin{1}=varargin{1}(temp);
varargin{2}=varargin{2}(temp);
varargin{3}=unique(varargin{3});

[varargout{:}] = interp1(varargin{:});

end

