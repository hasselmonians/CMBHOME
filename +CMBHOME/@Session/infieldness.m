function [ results ] = infieldness(varargin)
%infieldness Measurement of the "in-fieldness" - reflects how "in field"
%an animal is along its trajectory.
%
% ARGUMENTS
%   cel             n x 2 vector like [tetrode index, cell index].  If
%                   unassigned, uses self.cel. If has more than one cell,
%                   all returns are as cell arrays with the highest
%                   dimension as the cell.
%   params          (see below)
%
% RETURNS
%   results         A struct contaning the following fields
%                       flns    The in-fieldness at each time in
%                               infieldness_epoch.
%                       ts      The time points samples for flns.
%
% OPTIONAL PARAMETERS
%   fieldness_func  How flns is calculated. Can be any of the following
%                   options (defaults to 'Rate map'):
%                       Rate map    looks up the value from the rate map
%                                   with binside and gaussian smthing by
%                                   smth.  If this is chosen, results
%                                   will include
%                                   no_occupancy    as for RATEMAP
%                                   xdim            as for RATEMAP
%                                   ydim            as for RATEMAP
%                                   flns_map        Ratemap normalized by
%                                                   percentile of rate
%                                   rate_map         Ratemap before
%                                                   normailzation
%                                   flns_raw        flns from raw rate map
%                       Instantaneous
%                       spike rate  the the instantaneous spike rate with
%                                   bin of binside and gaussian smthing
%                                   by smth.
%                       Polar       looks up the value from the polar rate
%                                   map with bin binside and gaussian
%                                   smthing by smth If this is chosen,
%                                   results will include:
%                                   theta           as for
%                                                   DIRECTIONALTUNINGFCN
%                                   flns_map        Ratemap normalized by
%                                                   percentile of rate
%                                   rate_map         Ratemap before
%                                                   normailzation
%                                   flns_raw        flns from raw rate map
%   binside         The size of the bin for calculating flns.
%   smth          The amount of smthing for the calculation of
%                   flns.
%   speed_correction If 'arc', Resamples along the trajectory of the animal,
%                   changing sampling from a time frequency (s^-1) to a
%                   movememnt frequency (cm^-1).  If 'hd', resamples along
%                   the head tuning path of the animal, changing sampling
%                   from a time frequency to a turning frequency (deg^-1).
%                   Defaults to 'arc'. If this is chosed, results will
%                   contain:
%                       ts_raw          timestamps from raw flns
%                                       calculation (before speed
%                                       correction)
%                       cs              Distance animal ran
%  direction
%
%   The defaults for binside, smth, and speed_correction are dependant on
%   flns function.
%       --------------------------------------------------------
%       |                 |         |  Instantaneous |         |
%       |                 |Rate map |  Spike rate    | Polar   |
%       |-----------------+---------+----------------+---------|
%       |         binside |  1cm    |     0.5 sec    |  6 deg  |
%       |          smth |  5      |      0         |   1     |
%       |speed_correction | 'arc'   |    'none'      | 'polar' |
%       --------------------------------------------------------
%
%
%
%   infieldness_epoch The epoch used in calculating the flns. If
%                   not supplied, uses [0 inf].
%   subsample       The amount at which the movement data is subsampled for
%                   the calculation of flns. Defaults to 1 (no subsampling)
%   fc              How often along the trajectory it is subsampled
%                   during speed correction. Defaults to sampling as many
%                   times as is time sampled (if fc is empty).
%
% See also PASS_INDEX
%
% Jason R. Climer Revision: 1.1 Date 10/02/11

%% Parse Input

import CMBHOME.Utils.*;

p = inputParser;
p.StructExpand = true;
p.KeepUnmatched = true;

[self,cel,varargin]=infieldness_p(varargin{:});

p.addParamValue('fieldness_func','Rate map',@ischar);

p.parse(varargin{:});

fieldness_funcs = {
    'user'
    'rate map'
    'polar'
    'instantaneous spike rate'
    };

i = find(ismember(fieldness_funcs,lower(p.Results.fieldness_func)));

switch lower(p.Results.fieldness_func)
    case 'rate map'
        p.addParamValue('rate_map',[],@isnumeric);
        p.addParamValue('xdim',[],@isnumeric);
        p.addParamValue('ydim',[],@isnumeric);
        p.addParamValue('no_occupancy',[],@islogical);
    case 'polar'
        p.addParamValue('rate_map',[],@isnumeric);
        p.addParamValue('theta',[]);
end

if isempty(i), i=1;end;

binsides = {
    1
    1% cm
    3% degrees
    2 % seconds
    };

smths = {
    1
    5
    6
    0.0
    };

speed_corrections = {
    'none'
    'arc'
    'polar'
    'none'
    };

p.addParamValue('speed_correction',speed_corrections{i},@(x)any(strcmpi(x,speed_corrections)));
p.addParamValue('fc',[],@(x)isscalar(x)|isempty(x));

p.addParamValue('binside',binsides{i},@(x)isscalar(x)&x>0);
p.addParamValue('smth',smths{i},@(x)isscalar(x)&x>=0);

p.addParamValue('subsample',1,@(x)isscalar(x)&x>0&x==floor(x));
p.addParamValue('infieldness_epoch',[0 inf],@(x)size(x,2)==2);

p.addParamValue('dir',NaN);

p.parse(varargin{:});

% Get everything
subsample = p.Results.subsample;
binside = p.Results.binside;
smth = p.Results.smth;
dir = p.Results.dir;
infieldness_epoch = p.Results.infieldness_epoch;
fieldness_func = p.Results.fieldness_func;
speed_correction = p.Results.speed_correction;
fc = p.Results.fc;

if strcmpi(fieldness_func,'polar')&&~strcmpi(speed_correction,'none')&&isempty(fc)
    fc = 1/binside;
end

%% Setup values

results = struct;

% Set epochs

self.epoch = infieldness_epoch;
infieldness_epoch = self.epoch;

%% Handle multiple epochs
if size(infieldness_epoch,1)>1
    p = p.Results;
    results = struct;
    for i=1:size(fieldness_epcoh,1)
        p.infieldness_epoch = infieldness_epoch(i,:);
        k = infieldness(p);
        for j=fields(k)
            p.(j){i}=k(j);
        end
    end
else
    %% Calculate flns by fieldness_func
    
    % Assign cel
    self.cel = cel;
    
    if strcmpi(fieldness_func,'Rate map')
        
        % Get rate map
        [rate_map, xdim, ydim,~, no_occupancy] = self.sRateMap(cel,'continuize_epochs',1,'binside',binside,'std_smooth_kernel',smth);
        if ~isempty(p.Results.rate_map)
            rate_map = p.Results.rate_map;
        end
        if ~isempty(p.Results.xdim)
            xdim = p.Results.xdim;
        end
        if ~isempty(p.Results.ydim)
            ydim = p.Results.ydim;
        end
        if ~isempty(p.Results.no_occupancy)
            no_occupancy = p.Results.no_occupancy;
        end
        
        
        rate_map(no_occupancy) = NaN;
        
        rm_lookup = sort(rate_map(~isnan(rate_map)));
        rm_lookup = [rm_lookup linspace(0,1,numel(rm_lookup))'];
        
        [~,flns_map] = min(abs(bsxfun(@minus,rm_lookup(:,1),rate_map(~isnan(rate_map))')));
        temp = rate_map;
        temp(~isnan(temp)) = rm_lookup(flns_map(:),2);
        flns_map = temp;
        flns_map(no_occupancy) = NaN;
        
        % Find nearest bins
        if length(xdim)>1
            [~,i] = min(abs(bsxfun(@minus,xdim',self.sx(1:subsample:end)')));
            i = i';
        else
            i = ones(size(self.sx(1:subsample:end)));
        end
        if length(ydim)>1
            [~,j] = min(abs(bsxfun(@minus,ydim',self.sy(1:subsample:end)')));
            j = j';
        else
            j = ones(size(self.sy(1:subsample:end)))*ydim;
        end
        
        
        flns = flns_map(sub2ind(size(flns_map),j,i));
        flns_raw = rate_map(sub2ind(size(rate_map),j,i));
        
        results.no_occupancy = no_occupancy;
        results.xdim = xdim;
        results.ydim = ydim;
        results.flns_map = flns_map;
        results.rate_map = rate_map;
        results.flns_raw = flns_raw;
        
    elseif strcmpi(fieldness_func,'Polar')
        
        % Get polar rate map
        [rate_map, theta] = self.DirectionalTuningFcn(cel,'binsize',binside,'Continuize',1);
        
        rate_map = smooth(rate_map,smth/binside);
        rate_map(rate_map<0) = 0;
                
        if ~isempty(p.Results.rate_map)
            rate_map = p.Results.rate_map;
        end
        if ~isempty(p.Results.theta)
            theta = p.Results.theta;
        end
        
        rm_lookup = sort(rate_map(~isnan(rate_map)));
        rm_lookup = [rm_lookup linspace(0,1,numel(rm_lookup))'];
        
        [~,flns_map] = min(abs(bsxfun(@minus,rm_lookup(:,1),rate_map(:)')));
        flns_map = reshape(rm_lookup(flns_map,2),size(rate_map));
        
        % Find nearest bins
        hd = self.headdir(1:subsample:end);
        
        [~,tinds] = min(abs(bsxfun(@minus,theta,hd')));
        
        flns = flns_map(tinds);
        raw_flns = rate_map(tinds);
        
        results.rate_map = rate_map;
        results.flns_map = flns_map;
        results.theta = theta;
        results.raw_flns = raw_flns;
        
    elseif strcmpi(fieldness_func,'Instantaneous spike rate')
        
        maxl = 20000;
        i = 1;
        
        flns = [];
        inds = (1:subsample:length(self.ts))';
        
        while i<length(inds)
            top = i+maxl;
            if (top>length(inds)), top=length(inds);end;
            inds1 = inds(i:top);
            flns = [flns;sum(abs(repmat(self.spk.ts',[length(self.ts(inds1)) 1])-repmat(self.ts(inds1),[1 length(self.spk.ts')]))<=binside/2,2)/binside];
            i=top+1;
        end
        flns = CMBHOME.Utils.smthMat(flns,[0 5*smth/binside]);
    end
    
    if size(flns,1)==1
        flns = flns';
    end
    
    ts = self.ts(1:subsample:end);
    
    %% Interpolate over nans
    
    flns = nanInterp(flns);
    
    %% Speed Correction
    if ~strcmpi(speed_correction,'none')
        % Build arc
        results.ts_raw = ts;
        
        if strcmpi(speed_correction,'arc')
            if ~isnan(dir)
                vel = abs(diff(self.sx*cos(dir)+self.sy*sin(dir)));
                vel(isnan(vel))=0;
                arcs = cumsum([0;vel]);
                arcs = arcs(1:subsample:end);
            else
                vel = sqrt(diff(self.sx).^2+diff(self.sy).^2);
                vel(isnan(vel))=0;
                arcs = [0;cumsum(vel)];
                arcs = arcs(1:subsample:end);
            end
        elseif strcmpi(speed_correction,'polar')
            hds = self.headdir(1:subsample:end);
            hds(isnan(hds)) = 0;
            arcs = cumsum([0;abs(anglediff(hds,'r',180))]);
            arcs(isnan(hds)) = NaN;
            arcs = nanInterp(arcs);
        end
        arcs = arcs-arcs(1);
        
        % eliminate duplicate arcs
        [~,i] = unique(arcs,'first');
        [~,j] = unique(arcs,'last');
                
        dups = [i(~ismember(i,j)) j(~ismember(j,i))];
        arcs = arcs(i);
        
        removed = 0;
        
        for j=dups'
            i = j-removed;
            flns = [flns(1:i(1)-1); mean(flns(i(1):i(2))); flns(i(2)+1:end)];
            ts = [ts(1:i(1)-1); mean(ts(i(1):i(2))); ts(i(2)+1:end)];
            removed = removed+diff(i);
        end
        
        if isempty(fc)
            fc = length(flns)/max(arcs);
        end
        
        arc_axis = (0:1/fc:max(arcs))';
        while (arc_axis(end)<max(arcs))
            arc_axis = [arc_axis;arc_axis(end)+fc];
        end
        
        flns = interp1(arcs,flns,arc_axis,'nearest');
        ts = interp1(arcs,ts,arc_axis,'nearest');
        
        arc_axis = arc_axis(~isnan(flns));
        ts = ts(~isnan(flns));
        flns = flns(~isnan(flns));
        
        results.cs=arc_axis(1:length(flns));
    end
    
    results.ts = ts;
    results.flns = flns;
    
end
end

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

