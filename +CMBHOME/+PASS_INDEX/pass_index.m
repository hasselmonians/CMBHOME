function results = pass_index( varargin )
% PASS_INDEX - Calcultes the pass index and plots
%
% Calculates the pass index for the data passed. This function uses the
% pass_index_parser to generate its input structure.
%
% RESULTS = PASS_INDEX_PARSER(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG)
% RESULTS = PASS_INDEX_PARSER(POS_TS,POS,SPK_TS,LFP_TS,LFP_SIG,PARAMS)
%
%   ARGUMENTS
%   * POS_TS: Vector of time stamps for the sample state
%   * POS: MXN matrix of the sample state, where M is the number of samples
%   and N is the dimensions of POS
%   * SPK_TS: Spike times for the cell
%   * LFP_TS: Time stamps for the local field potential (LFP) Sample
%   * LFP_SIG: The LFP signal
%
%   PARAMETERS
%   * plots: Default false. If false, doesn't plot. If true or 'all', plots
%   all possible plots. Can be true,'all',any or or a cell array of the
%   following:
%     	�Trajectory�: (1-3D only) This plots the trajectory of the animal with the spikes of the cell colored by the pass index of the spike.
%       �Rate map�: (1-3D only) This plots the rate map of the cell.
%   	�Field index map�: (1-3D only) This plots the field index map of the cell. See table 1 for notes on custom implementation.
%   	�Scatter plot�: Shows the pass index versus two cycles of the lfp phase for all spikes, and calculates the linear circular correlation using the circular-linear correlation19. This can be calculated using the [corr_val,p,s,b]=kempter_lincirc(linear,circular) function included, and are output in the results struct (see 5.4 and Table 1).
%   	�Density map�: Shows the density of the scatter plot, using 100 phase bins and 40 pass index bins, then smoothing with a pseudo-Gaussian kernel with width of 1.5 pixel.
%   * subplots: Default []. If the parameter �subplots� is set to a nX2
%   matrix, where n is the number of plots, each of the rows of the
%   subplots parameter will  be used to determine where the subplot will
%   go. Otherwise, will plot as squarly as possible, biased wider than
%   taller.
%
%   PASS INDEX PARAMETERS
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
%   * RESULTS: A struct containing the following fields:
%       rate_map	Occupancy normalized rate map for the cell.
%       field_index_map	Field index map as calculated by the field index function.
%       centers the centers for the maps.
%       ts	Resampled time stamps.
%       cs	Resampled axis positions. As a default, the distance along the trajectory (cm).
%       field_index	Field index at every point along the resampled axis.
%       filtered_field_index Field index after filter is applied
%       pass_index	Pass index at every point along the resampled axis.
%       spk_pass_index	The pass index at the time of each spike.
%       spk_theta_phase	The theta phase at the time of each spike.
%       rho	The linear-circular correlation coefficient from kempter_lincirc.
%       p	The significance level of the correlation.
%       s	The slope in � per pass.
%       is_precessing	True if p<0.05 and -1440<s<-22 � per pass.
%       density	Density map, with LFP phase broken into 100 bins and pass index broken into 40.
%       filtered_lfp	The filtered LFP
%       filtered_lfp_phase	The phase of the filtered LFP
%
%  From pass_index. Release 2013-09-13 v0.1 from Jason Climer
%  (jason.r.climer@gmail.com)
% Parse inputs
import CMBHOME.PASS_INDEX.*;
import CMBHOME.Utils.*;

P = pass_index_parser(varargin{:});
for i = fields(P.Results)'
    eval([i{1} ' = P.Results.' i{1} ';']);
end

% Calculate field index
if isequal(class(field_index),'function_handle')
    fi_fun = field_index;
    try
        field_index = fi_fun(varargin{:});
    catch
        error('Pass_index failed. Did you FixPos?')
    end
    if any(cellfun(@(x)isequal(x,'field_index'),varargin))
        varargin{find(cellfun(@(x)isequal(x,'field_index'),varargin))+1}=field_index;
    else
        varargin = [varargin {'field_index',field_index}];
    end
end

% Resample
if isequal(class(sample_along),'function_handle')
    sample_along = sample_along(varargin{:});
    if ~any(cellfun(@(x)isequal(x,'sample_along'),varargin))
        varargin = [varargin {'sample_along',sample_along}];
    else
        varargin{find(cellfun(@(x)isequal(x,'sample_along'),varargin))+1}=sample_along;
    end
end

% Filter
[filtered] = filter_band(varargin{:});

% Calculate Pass Index
pass_index = angle(hilbert(filtered))/pi;
spk_pass_index = (mod(interp1(sample_along(:,2),unwrap(pass_index*pi),spk_ts,'nearest','extrap')+pi,2*pi)-pi)/pi;

% Calculate LFP
if ~isempty(lfp_sig)
    [filtered_lfp,lfp_phase] = lfp_filter(varargin{:});
    spk_theta_phase = mod(interp1(lfp_ts,unwrap(lfp_phase),spk_ts)+pi,2*pi)-pi;
end

%% Packing results
results = struct();

[map,centers,occupancy] = rate_map(varargin{:});
results.rate_map = map;
results.centers = centers;
results.occupancy = occupancy;

fi_map = NaN;

try
    fi_map = fi_fun(varargin{:},'get_map',true);
    if ~isequal(size(fi_map),size(map))
        fi_map = NaN;
    end
catch err
end

results.field_index_map = fi_map;
results.ts = sample_along(:,2);
results.cs = sample_along(:,1);
results.field_index = sample_along(:,3);
results.pass_index = pass_index;
results.filtered_field_index = filtered;
results.spk_pass_index = spk_pass_index;
results.spk_theta_phase = spk_theta_phase;
results.filtered_lfp = filtered_lfp;
results.filtered_lfp_phase = lfp_phase;

[rho,p,s,b] = kempter_lincirc(spk_pass_index,spk_theta_phase); % Correlation
results.rho = rho;
results.p = p;
results.s = s;
results.b = b;
results.is_precessing = p<0.05&&rad2deg(pi*s)<-22&&rad2deg(pi*2)>-1440;

% Density
dens_oc = histcn([interp1(pos_ts,pass_index,lfp_ts,'nearest')' mod(lfp_phase,2*pi)'],linspace(-1,1,41),linspace(0,2*pi,101));
dens_oc = dens_oc*mean(diff(lfp_ts));
[density,~,dens_centers] = histcn([spk_pass_index mod(spk_theta_phase,2*pi)],linspace(-1,1,41),linspace(0,2*pi,101));
density = density./dens_oc;
h = 1.5;
myfilter = fspecial('gaussian',[4 4]*h, h);
density = imfilter(density,myfilter,'replicate');
results.density = density;
results.dens_centers = dens_centers;

%% Plotting
% Parse plot stuff
ip = plot_pass_index_parser(varargin{:},'results',results);
%% Plot
if ~PLOT_OVERRIDE&&numel(ip.Results.plots)>0
    plot_pass_index(varargin{:},'results',results);
end
end

