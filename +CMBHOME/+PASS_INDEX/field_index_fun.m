function [ field_index ] = field_index_fun( varargin )
% FIELD_INDEX_FUN - Calculated the field index.
%
% Calculates the field index along a trajectory. This function uses the
% pass_index_parser to generate its input structure.
%
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS)
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG)
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG,PARAMS)
% FIELD_INDEX = FIELD_INDEX_FUN(POS_TS,POS,SPK_TS,[],[],PARAMS)
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
%   * get_map: Default false. If true, returns a map instead of sampling
%   along the trajectory
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
%   * FIELD_INDEX: If method is 'grid', then the bins of the rate map are 
%   percentile ranked between 0 and 1. Thus, a field index value of 0.5 
%   indicates that 50% of the bins have a lower rate that the bin in 
%   question. If method is 'place', then teh bins of the rate map are 
%   linearly normalized between 0 and 1.
%
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com) 
import CMBHOME.PASS_INDEX.*;

p = pass_index_parser(varargin{:});
for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end

p.addParamValue('get_map',false);
p.parse(varargin{:});
get_map = p.Results.get_map;

[ map,centers,occupancy,~,~ ] = rate_map( varargin{:} );

fi_map = map;
fi_map(occupancy==0) = NaN;
switch method
    case 'grid'
        [fi_dist,i] = sort(fi_map(:));
        if any(isnan(fi_dist(:)))
            k = 1:find(isnan(fi_dist),1);
        else
            k = 1:numel(fi_dist);
        end
        fi_dist(k) = linspace(0,1,numel(k));
        fi_map(i) = fi_dist;
    case 'place'
        fi_map = (fi_map-nanmin(fi_map(:)))/range(fi_map(:));
end

if (get_map)
    field_index = fi_map;
else

    pos2 = cellfun(@(x,y)round((x-min(y))/binside),cellfun(@(x)pos(:,x),num2cell(1:size(pos,2)),'UniformOutput',false),centers,'UniformOutput',false);
    field_index = zeros(size(pos,1),1);
    field_index(any([pos2{:}]<=0|[pos2{:}]>repmat(size(map),[numel(pos_ts) 1]),2)) = NaN;
    pos2 = cellfun(@(x)x(~isnan(field_index)),pos2,'UniformOutput',false);
    field_index(~isnan(field_index)) = fi_map(sub2ind(size(map),pos2{:}));
    field_index = interp1(pos_ts(~isnan(field_index)),field_index(~isnan(field_index)),pos_ts,'linear','extrap');
end 
   

