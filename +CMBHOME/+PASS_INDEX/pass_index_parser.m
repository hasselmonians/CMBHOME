function p = pass_index_parser(varargin)
% PASS_INDEX_PARSER - Parses pass-index related inputs
%
% Parses pass index related inputs for pass_index and function that 
% pass_index uses
%
% P = PASS_INDEX_PARSER(POS_TS,POS,SPK_TS)
% P = PASS_INDEX_PARSER(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG)
% P = PASS_INDEX_PARSER(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG,PARAMS)
%
%   ARGUMENTS
%   * POS_TS: Vector of time stamps for the sample state
%   * POS: MXN matrix of the sample state, where M is the number of samples
%   and N is the dimensions of POS
%   * SPK_TS: Spike times for the cell
%
%   OPTIONAL ARGUMENTS
%   * LFP_TS: Time stamps for the local field potential (LFP) Sample
%   * LFP_SIG: The LFP signal
%
%   PARAMETERS
%   * method: Default 'grid'. Can be 'grid','place', or custom. Updates
%   other unset fields for these techniques.
%   * binside: Default 2*N, where N is the dimensionality of POS. Side of
%   the bins for rate mapping.
%   * smth_width: Default 3*BINSIDE, width of Gaussian smoothing kernel
%   * field_index: Default @field_index_fun, can be a vector of the same 
%   number of elements as pos_ts, or can be a function handle which takes 
%   in the same parameters as pass_index.
%   * sample_along: Default 'auto', can be 'arc_length', 'raw_ts', or a 
%   nX2 matrix where n is the number of resampled steps, the first column 
%   is the resampled timestamps and the second column is the sampled values, 
%   or a function handle that returns a nX2 matrix as described above. Set
%   from 'auto' to 'arc_length' if method is 'place' or 'grid'.
%   * filter_band: Default 'auto', can be any positive frequency range in 
%   cycles/unit sampled along using the �filter_band� parameter. 
%   Additionally, filter_band can be a function handle which returns a 
%   modified signal. Set from 'auto' to [0.0749 0.0029] if 'method' is
%   'grid' and to the [3*D 1/6*D].^-1, where D is the field width 
%   determined by finding the N-dimensional volume of the region with at 
%   least 10% of the maximum firing rate, and calculating the diameter of 
%   the n-ball with the same volume. 
%   * lfp_filter: Default [6 10]. can be changed to any frequency range in 
%   Hz as [low high] or as a function handle with the form lfp_phases = 
%   custom_phase_func(lfp_ts,lfp_sig) for custom phase estimation, for 
%   example, by taking asymmetry into account
%   * lfp_subsample_fs: Defaults 500.
%
%   RETURNS
%   * P: An input parser with all inputs correctly parsed and default
%   values updated.
%
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com) 

%% Build Parser
import CMBHOME.PASS_INDEX.*;

p = inputParser;
p.KeepUnmatched = true;

% Get data
p.addParameter('pos_ts',@(x)isnumeric(x)&sum(size(x)==1)>=1);
p.addParameter('pos',@isnumeric);
p.addParameter('spk_ts',@(x)isnumeric(x)&sum(size(x)==1)>=1);
p.addParameter('lfp_ts',[],@(x)isempty(x)||(isnumeric(x)&&sum(size(x)==1)>=1));
p.addParameter('lfp_sig',[],@(x)isempty(x)||(isnumeric(x)&&sum(size(x)==1)>=1));
p.parse(varargin{:});
for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end

% interp over NaNs
for i=1:size(pos,2)
   pos(:,i) = interp1(pos_ts(~isnan(pos(:,i))),pos(~isnan(pos(:,i)),i),pos_ts,'linear','extrap');
end
varargin{4} = pos;

% Add parser terms
p.addParamValue('method','grid',@(x)(ischar(x)&&...
    ismember(x,{'custom','grid','place'}))...
    ||isequal(class(x),'function_handle'));
p.addParamValue('binside','auto',@(x)isscalar(x)||isequal(x,'auto'));
p.addParamValue('smth_width','auto',@(x)isscalar(x)||isequal(x,'auto'));
p.addParamValue('field_index',@field_index_fun,...
    @(x)isequal(class(x),'function_handle')||...
    (isvector(x)&&isnumeric(x)&&numel(x)==numel(pos_ts)));
p.addParamValue('sample_along','auto',@(x)(ischar(x)&&...
    ismember(x,{'auto','arc_length','raw_ts'}))||...
    (isnumeric(x)&&size(x,2)==3&&all(diff(x(:,1))>0))||...
    isequal(class(x),'function_handle'));
p.addParamValue('filter_band','auto',...
    @(x)(ischar(x)&&ismember(x,{'auto'}))||...
    (isnumeric(x)&&numel(x)==2&&x(1)<x(2))||...
    isequal(class(x),'function_handle'));
p.addParamValue('lfp_filter',[6 10],@(x)(isnumeric(x)&&isequal(size(x),[2 1]))||...
    isequal(class(x),'function_handle'));
p.addParamValue('lfp_subsample',500,@(x)isnumeric(x)&&isscalar(x));
p.addParamValue('PLOT_OVERRIDE',0);
p.parse(varargin{:});

%%
% lfp_subsample
lfp_ts = p.Results.lfp_ts(~isnan(p.Results.lfp_sig));
lfp_sig = p.Results.lfp_sig(~isnan(p.Results.lfp_sig));

if numel(unique(diff(lfp_ts)))>1||unique(diff(lfp_ts))^-1~=p.Results.lfp_subsample
    temp = lfp_ts;
    lfp_ts = lfp_ts(1):1/p.Results.lfp_subsample:lfp_ts(end);
    lfp_sig = interp1(temp,lfp_sig,lfp_ts,'linear','extrap');
end

% Replace in VARARGIN
varargin{8} = lfp_ts;
varargin{10} = lfp_sig;
p.parse(varargin{:});

%% Update auto scalars
if isequal(p.Results.binside,'auto')
    switch p.Results.method
        case 'grid'
            if ismember('binside',p.UsingDefaults)
                varargin = [varargin {'binside',2*size(pos,2)}];
            else
                varargin{find(cellfun(@(x)isequal(x,'binside'),varargin))+1}=2*size(pos,2);
            end
        case 'place'
            if ismember('binside',p.UsingDefaults)
                varargin = [varargin {'binside',2*size(pos,2)}];
            else
                varargin{find(cellfun(@(x)isequal(x,'binside'),varargin))+1}=2*size(pos,2);
            end
    end
    p.parse(varargin{:});
end

if isequal(p.Results.smth_width,'auto')
    switch p.Results.method
        case 'grid'
            if ismember('smth_width',p.UsingDefaults)
                varargin = [varargin {'smth_width',3*p.Results.binside}];
            else
                varargin{find(cellfun(@(x)isequal(x,'smth_width'),varargin))+1}=3*p.Results.binside;
            end
        case 'place'
            if ismember('smth_width',p.UsingDefaults)
                varargin = [varargin {'smth_width',3*p.Results.binside}];
            else
                varargin{find(cellfun(@(x)isequal(x,'smth_width'),varargin))+1}=3*p.Results.binside;
            end
    end
    p.parse(varargin{:});
end

%% Update filter band by the method
if isequal(p.Results.filter_band,'auto')
    if ischar(p.Results.method)
        switch p.Results.method
            case 'grid'
                
                if ismember('filter_band',p.UsingDefaults)
                    varargin = [varargin {'filter_band',[2*170 26.7/2].^-1}];
                else
                    varargin{find(cellfun(@(x)isequal(x,'filter_band'),varargin))+1}=[2*170 26.7/8].^-1;
                end
                
            case 'place'
                map = rate_map(varargin{1:5},'binside',p.Results.binside,'smth_width',p.Results.smth_width);
                V = sum(map(:)>0.1*max(map(:)))*p.Results.binside^size(pos,2);
                n = size(pos,2);
                k = floor(n/2);
                
                if mod(n,2) % odd
                    r = (prod(1:2:n)*V/(2^(k+1)*pi^k))^(1/(2*k+1));
                else % even
                    r = (factorial(k)*V)^(1/(2*k))/sqrt(pi);
                end
  
                if ismember('filter_band',p.UsingDefaults)
                    varargin = [varargin {'filter_band',[r*6 r/3].^-1}];
                else
                    varargin{find(cellfun(@(x)isequal(x,'filter_band'),varargin))+1}=[r*6 r/3].^-1;
                end
        end
    end
    p.parse(varargin{:});
end

% Replace with function handle
if isnumeric(p.Results.filter_band)&&numel(p.Results.filter_band)==2&&p.Results.filter_band(1)<p.Results.filter_band(2)
    band = p.Results.filter_band;
    varargin{find(cellfun(@(x)isequal(x,'filter_band'),varargin))+1}=...
        @(varargin)if_filter_func(band,varargin{:});
    p.parse(varargin{:});
end

% Replace with function handle
if isnumeric(p.Results.lfp_filter)&&numel(p.Results.lfp_filter)==2&&p.Results.lfp_filter(1)<p.Results.lfp_filter(2)
    band = p.Results.lfp_filter;
    if ismember('lfp_filter',p.UsingDefaults)
        varargin = [varargin {'lfp_filter',@(varargin)lfp_filter_fun(band,varargin{:})}];
    else
        varargin{find(cellfun(@(x)isequal(x,'lfp_filter'),varargin))+1}=...
            @(varargin)lfp_filter_fun(band,varargin{:});
    end
    p.parse(varargin{:});
end

%% Update sample along by method
if isequal(p.Results.sample_along,'auto')
    if ischar(p.Results.method)
        switch p.Results.method
            case {'grid','place'}
                
                if ismember('sample_along',p.UsingDefaults)
                    varargin = [varargin {'sample_along','arc_length'}];
                else
                    varargin{find(cellfun(@(x)isequal(x,'sample_along'),varargin))+1}=[2*170 26.7/8].^-1;
                end
        end
    end
    p.parse(varargin{:});
end

% Replace chars with function handles
if ischar(p.Results.sample_along)
    switch p.Results.sample_along
        case 'arc_length'
            if ismember('sample_along',p.UsingDefaults)
                varargin = [varargin {'sample_along',@sample_along_arc}];
            else
                varargin{find(cellfun(@(x)isequal(x,'sample_along'),varargin))+1}=@sample_along_arc;
            end
        case 'raw_ts'
            if ismember('sample_along',p.UsingDefaults)
                varargin = [varargin {'sample_along',@sample_along_ts}];
            else
                varargin{find(cellfun(@(x)isequal(x,'sample_along'),varargin))+1}=@sample_along_ts;
            end
    end
    p.parse(varargin{:});
end


end

function [filtered,phase] = lfp_filter_fun(band,varargin) % Filters LFP
import CMBHOME.PASS_INDEX.*;
p = pass_index_parser(varargin{:});
for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end
Fs = mode(diff(lfp_ts)).^-1;
Wn = band/(Fs/2);
[b,a] = butter(3,Wn);
filtered = filtfilt(b,a,lfp_sig);
phase = angle(hilbert(filtered));
end

function [filtered] = if_filter_func(band,varargin) % Basic field index filter - resamples and filteres
import CMBHOME.PASS_INDEX.*;
p = pass_index_parser(varargin{:});
for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end

if isequal(class(p.Results.sample_along),'function_handle')
    sample_along = p.Results.sample_along(varargin{:});
end

cc = sample_along(:,1);
ts2 = sample_along(:,2);
resampled = sample_along(:,3);

Fs = mean(diff(cc)).^-1;
Wn = band/(Fs/2);
Wn(Wn>=1) = 1-1*10^-5;
Wn(Wn<0) = 0;
[b,a] = butter(3,Wn);
filtered = filtfilt(b,a,resampled);

end

function [sample_along] = sample_along_arc(varargin) % Samples along arc traversed
import CMBHOME.PASS_INDEX.*;
p = pass_index_parser(varargin{:});
for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end

arc = cumsum([0;sqrt(nansum((diff(pos)).^2,2))]);
ts2 = pos_ts(diff(arc)>0);
arc = arc(diff(arc)>0);
cc = linspace(0,max(arc),size(pos,1))';
ts2 = interp1(arc,ts2,cc);
resampled = interp1(pos_ts,field_index,ts2);

sample_along = [cc,ts2,resampled];

end

function [sample_along] =  sample_along_ts(varargin) % Passes everything back as it was
import CMBHOME.PASS_INDEX.*;
p = pass_index_parser(varargin{:});
for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end

ts2 = pos_ts;
cc = pos_ts;
flns2 = field_index;

sample_along = [cc,ts2,resampled];

end