function [ self, cel, args] = infieldness_p(varargin)
% Handles infieldness required fields to an input parser

if strcmp(class(varargin{1}),'inputParser')
    p = varargin{1};
    varargin = varargin(2:end);
else
    p = inputParser;
end
p.StructExpand = true;
p.KeepUnmatched = true;

if isstruct(varargin{1})
   if isfield(varargin{1},'cel')
       [self, cel, args] = infieldness_p(varargin{1}.self,varargin{1}.cel,varargin{:});
   else
       [self, cel, args] = infieldness_p(varargin{1}.self,varargin{:});
   end
else
    p.addRequired('self',@(x)isstruct(x)||strcmp(class(x),'CMBHOME.Session'));
    p.parse(varargin{1});
    self = p.Results.self;
    p.addOptional('cel',self.cel,@(x)isempty(x)||all(ismember(x,self.cells,'rows')));
    p.parse(varargin{:});
    self = p.Results.self;
    cel = p.Results.cel;
    args = varargin(find(cellfun(@ischar,varargin), 1 ):end);
end

end

