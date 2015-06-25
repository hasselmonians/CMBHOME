function [ map,centers,occupancy,pos_loc,spk_loc ] = rate_map( varargin )
% RATE_MAP - Calculated the occupancy normalized rate map
%
% Calculates the N-dimensional occupancy normalized rate map. rate_amp 
% uses the pass_index_parser inputs.
%
% [MAP,CENTERS,OCCUPANCY,POS_LOC,SPK_LOC] = RATE_MAP(POS_TS,POS,SPK_TS)
% [MAP,CENTERS,OCCUPANCY,POS_LOC,SPK_LOC] = RATE_MAP(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG)
% [MAP,CENTERS,OCCUPANCY,POS_LOC,SPK_LOC] = RATE_MAP(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG,PARAMS)
% [MAP,CENTERS,OCCUPANCY,POS_LOC,SPK_LOC] = RATE_MAP(POS_TS,POS,SPK_TS,[],[],PARAMS)
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
%   * spkpos: Default []. If left empty, uses spk_pos to find spike
%   positions.
%
%   PASS_INDEX_PARSER PARAMETERS
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
%   cycles/unit sampled along using the ‘filter_band’ parameter. 
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
%
%   RETURNS
%   * MAP - occupancy normalized rate map
%   * CENTERS - cell array with the centers of the bins along each axis
%   * OCCUPANCY - time spent in each bin
%   * POS_LOC - index in map of each position tracked
%   * SPK_LOC - index in map of each spike tracked
%
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com) 
import CMBHOME.PASS_INDEX.*;
import CMBHOME.Utils.*;
% Parse inputs
p = pass_index_parser(varargin{:});
for i = fields(p.Results)'
eval([i{1} ' = p.Results.' i{1} ';']);
end

% Interpolate over missing/unevenly spaced data
pos_orig = pos;
pos_orig_ts = pos_ts;
temp = pos_ts(~any(isnan(pos),2));
pos = pos(~any(isnan(pos),2),:);
pos_ts = linspace(min(temp),max(temp),numel(temp));
pos = interp1(pos_orig_ts,pos_orig,pos_ts,'linear','extrap');

% Spike positions
p2 = p;
p2.addParamValue('spkpos',[],@isnumeric);
p2.parse(varargin{:});
spkpos = p2.Results.spkpos;
if isempty(spkpos)&&~isempty(spk_ts)
    spkpos = spk_pos(pos_ts,pos,spk_ts);
end
  
% Build dimension cell array
dims = minmax(pos')';
dims(1,:) = floor(dims(1,:)/binside)*binside;
dims(2,:) = ceil(dims(2,:)/binside)*binside;
dims = cellfun(@(x,y)x:binside:y,num2cell(dims(1,:)),num2cell(dims(2,:)),'UniformOutput',false)';

% Find occupancy
[occupancy,~,centers] = histcn(pos,dims{:});
occupancy = occupancy*mean(diff(pos_ts));

% Count spikes
[map,~,~,spk_loc] = histcn(spkpos,dims{:});
map = map./occupancy; % occupancy normalize
map(occupancy==0) = NaN;

[~,~,~,pos_loc] = histcn(pos_orig,dims{:});

% % SMOOTHING % %

% Patch inner occupancy holes with the mean rate of the surrounding pixels
props = regionprops(bwlabeln(isnan(map)),'PixelIdxList','PixelList');
temp = cellfun(@(x)[',props(i).PixelList(:,' num2str(x) ')'],num2cell(1:size(pos,2)),'UniformOutput',false);
temp = [temp{:}];
temp = temp(2:end);
for i=1:numel(props)
    eval(['[' temp ']=ind2sub(size(map),props(i).PixelIdxList);']);
end
props = props(~cellfun(@(x)any(x(:)==1)||...
    eval([repmat('any(',[1 size(pos,2)]) 'x==repmat(size(map),[size(x,1) 1])' repmat(')',[1 size(pos,2)])]),...
    {props.PixelList}));
temp = cellfun(@(x)[',y(:,' num2str(x) ')'],num2cell(1:size(pos,2)),'UniformOutput',false);
temp = [temp{:}];
for i=1:numel(props)
   x = props(i).PixelList;
   pixelselect = false(size(occupancy));
   for j=1:size(x,2)
       y = x;
       y(:,j) = y(:,j)+1;
       eval(['pixelselect(sub2ind(size(occupancy)' temp '))=true;']);
       y(:,j) = y(:,j)-2;
       eval(['pixelselect(sub2ind(size(occupancy)' temp '))=true;']);       
   end
   y = x;
   eval(['pixelselect(sub2ind(size(occupancy)' temp '))=false;']);;
   map(props(i).PixelIdxList) = nanmean(map(pixelselect));
end

%% Make gaussian kernel & smooth
sigma = smth_width/binside/2;
M = floor(3*sigma);
m = repmat({-M:M},[size(pos,2) 1]);
h = cell(size(m));
[h{:}] = ndgrid(m{:});
k = mvnpdf(cell2mat(cellfun(@(x)x(:),h','UniformOutput',false)),zeros(size(pos,2),1)',eye(size(pos,2))*sigma);
if size(pos,2)>1
    k = reshape(k,repmat(M*2+1,[1 size(pos,2)]));
end
map(isnan(map)) = 0;
map = convn(map,k,'same');
map(occupancy==0) = 0;

end

