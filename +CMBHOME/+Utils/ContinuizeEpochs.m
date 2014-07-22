function varargout = ContinuizeEpochs(varargin)
% varargout = CMBHOME.Utils.ContinuizeEpochs(varargin);
%
% Converts all cell arrays passed through varagin to arrays in varargout
%
% If the cell array is Mx1, then that vararg output is a column vector. If
% the cell array is MxN, then the vararg output is a matrix max(cat(M) x N
% padded by zeros.
%
% This is useful in the CMBHOME framework: 

% takes all data from cell arrays in vargin, and continuizes them into
% vectors in vargout

varargout = cell(length(varargin),1);

for i = 1:length(varargin)
    
    if iscell(varargin{i}) && size(varargin{i},2)==1 % if only multiple epochs, but one cell
        varargout{i} = vertcat(varargin{i}{:});
    elseif iscell(varargin{i}) && size(varargin{i},2)>1 % if multiple epochs and cells
        
        cat_ca = cell(1, size(varargin{i},2));
        
        for ii = 1:size(varargin{i},2)
        
            cat_ca{1, ii} = vertcat(varargin{i}{:,ii});
            
        end
        
        varargout{i} = CMBHOME.Utils.CatCA(cat_ca);
        
    else
        varargout{i} = varargin{i};
    end
    
end