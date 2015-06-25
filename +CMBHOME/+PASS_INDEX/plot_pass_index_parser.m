function ip = plot_pass_index_parser( varargin )

import CMBHOME.PASS_INDEX.*;
if any(cellfun(@ischar,varargin))
    varargin1 = varargin(1:find(cellfun(@ischar,varargin),1)-1);
    varargin = varargin(find(cellfun(@ischar,varargin),1):end);    
else
    varargin1 = varargin;
    varargin = {};
end

all_plots = {'Trajectory','Rate map','Field index map','Scatter plot','Density map'};
ip = inputParser;
ip.KeepUnmatched=true;
ip.addParamValue('plots',true,@(x)(iscellstr(x)&&...
    all(ismember(x,all_plots)))||isequal(x,'all')||isequal(x,0)||isequal(x,1)||...
    (ischar(x)&&ismember(x,all_plots)));
ip.addParamValue('results',struct([]),@isstruct);
ip.parse(varargin{:});

if isempty(ip.Results.results)
    results = pass_index(varargin1{:},varargin{:},'PLOT_OVERRIDE',1);
    if ismember('results',ip.UsingDefaults)
        varargin = [varargin {'results',results}];
    else
        varargin{find(cellfun(@(x)isequal(x,'results'),varargin))+1}=...
            results;
    end
end

plots = ip.Results.plots;
if ~iscell(plots)
    switch plots
        case {1,'all'}
            plots = all_plots;
        case 0
            plots = {};
        otherwise
            plots = {plots};
    end
end
if ismember('plots',ip.UsingDefaults)
    varargin = [varargin {'plots',plots}];
else
    varargin{find(cellfun(@(x)isequal(x,'plots'),varargin))+1}=plots;
end

% keyboard

subplots = floor(sqrt(numel(plots)));
subplots = [subplots ceil(numel(plots)/subplots)];
ip.addParamValue('subplots',subplots);
ip.parse(varargin{:});
subplots = ip.Results.subplots;
if ~iscell(ip.Results.subplots)
   subplots = cellfun(@(x){ip.Results.subplots(1) ip.Results.subplots(2) x},num2cell(1:numel(plots)),'UniformOutput',false)'; 
   subplots = cat(1,subplots{:});
end
if ismember('subplots',ip.UsingDefaults)
    varargin = [varargin {'subplots',subplots}];
else
    varargin{find(cellfun(@(x)isequal(x,'subplots'),varargin))+1}=subplots;
end

ip.addParamValue('units','cm',@ischar);
ip.parse(varargin{:});

end

