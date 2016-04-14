function [ results ] = plot_phase_precession( self, varargin )
%PLOT_PHASE_PRECESSION Makes the phase precession plot
%
% ARGUMENTS
%
% OPTIONAL ARGUMENTS
%   cel     1 x 2 vector like [tetrode index, cell index].  If
%                   unassigned, uses self.cel.
%   lfp     Scalar integer that coresponds to which EEG channel is used. If
%           unassigned, uses self.active_lfp
%
% OPTIONAL PAREMETERS
%   plots   Cell array of which figures to plot.  By defuault, plots all,
%           but can be any of the following in any order.
%               Rate map    Plots the rate map.  Uses this instead of
%                           trajectory if fieldness_func is 'Polar'
%               Infieldness Plots the infieldness map.
%                   map
%               Trajectory  Plots the trajectory of the animal over the
%                           whole trial in grey, the selected magnitude
%                           ranges in black, and the spikes shown as the
%                           centered fieldness of each spike. Can be
%                           followed by..
%                           :Infieldness, default. Plots the infieldness of
%                           each spike
%                           :Rate Map.  Plots the rate map value of each
%                           spike
%                           :Pass Index Plots the pass index of each spike
%               Scatter     The scatter plot of three cell2mat(spk_theta) cycles over one
%                           pass through the field
%               Density     Plots the density plot of three cell2mat(spk_theta) cycles
%                           versus one pass through the field.
%               Pass Index Histogram   Plots a histogram with 20 bins of the position
%                           in field of the spikes.
%               Theta Histogram  Plots a histogram with 20 bins of the
%                           cell2mat(spk_theta) phases of the spikes
%               Correl  Plots the spike time autocorrellogram
%   Parameters are those of INFIELDNESS and PASS_INDEX. See
%   the help files for these functions for more information.
%
% RETURNS
%   results     Results of the underlying INFIELDNESS and PASS_INDEX
%
% See also FIELDNESS, FIELD_PASS_INDEX, SPK_FUNC
%
% Jason R. Climer Revision: 0.1 Date 08/08/11

if ~all(cellfun(@(x)isempty(x)||x.theta_phaseAppended,num2cell(self.b_lfp)))
   self.b_lfp = self.b_lfp.AppendTheta; 
end

import CMBHOME.Utils.*;

VALID_PLOTS = {'Auto','Scatter','Density','Histogram:Pass Index','Histogram:theta','Correl'};
OPTIONAL_PLOTS = {'Trajectory','Infieldness map','Infieldness map:polar','Infieldness map:spatial','Rate map:polar','Rate map','Rate map:spatial','Trajectory:Rate Map','Trajectory:Pass Index','Trajectory:Infieldness'};

p = inputParser;
p.addOptional('lfp',self.active_lfp,@(x)ismember(x,1:length(self.b_lfp)));
p.StructExpand = true;
p.KeepUnmatched = true;

[ self, cel, varargin ] = CMBHOME.Utils.infieldness_p(p,self,varargin{:});

p.addParamValue('plots',VALID_PLOTS,@(x)(ischar(x)&&(ismember(lower(x),cellfun(@lower,VALID_PLOTS,'UniformOutput',false))||ismember(lower(x),cellfun(@lower,OPTIONAL_PLOTS,'UniformOutput',false))))||...
    all(ismember(cellfun(@lower,x,'UniformOutput',false),cellfun(@lower,VALID_PLOTS,'UniformOutput',false))|ismember(cellfun(@lower,x,'UniformOutput',false),cellfun(@lower,OPTIONAL_PLOTS,'UniformOutput',false))));
p.addParamValue('traj_epoch',self.epoch,@(x)isnumeric(x)&size(x,2)==2);
p.addParamValue('fieldness_func','Rate map',@ischar);
p.addParamValue('infieldness_epoch',self.epoch);
p.addParamValue('pad',0.2);
p.addParamValue('siz','auto');
p.addParamValue('plot_inds','auto');
p.addParamValue('results',[]);


p.parse(self,varargin{:});

plots = p.Results.plots;
lfp = p.Results.lfp;
fieldness_func = p.Results.fieldness_func;
pad = p.Results.pad;
pad = [-pad pad];
traj_epoch = p.Results.traj_epoch;
infieldness_epoch = p.Results.infieldness_epoch;
siz = p.Results.siz;
plot_inds = p.Results.plot_inds;

if ~iscell(plots), plots = {plots}; end;
plots = cellfun(@lower,plots,'UniformOutput',false);

auto = '';

switch lower(fieldness_func)
    case 'rate map'
        auto = 'Trajectory:pass index';
    case 'polar'
        auto = 'Infieldness map:polar';        
    case 'instantaneous spike rate'
        auto = 'Trajectory:pass index';
end

if any(ismember(plots,'auto'))
    plots{ismember(plots,'auto')}=auto;
end

self.active_lfp = lfp;

if isempty(p.Results.results)
    results = pass_index(self,'cel',cel,varargin{:});
else
    results = p.Results.results;
end

% keyboard

try
    if isempty(self.b_lfp(self.active_lfp).theta_phase)
        throw(MException('a:b','a'));
    end
catch err
    throw(MException('plot_phase_precession:thetaErr','Theta not correctly imported'));
end

self.cel = cel;

self.epoch = [0 inf];

results.theta_cycle = [1;diff(self.lfp.theta_phase)<-pi];
results.theta_cycle = cumsum(results.theta_cycle);
results.lfpts = self.lfp.ts;

self.epoch = infieldness_epoch;

results.b_raw = zeros(size(self.b_ts));

k = self.ind;
if iscell(k), k=cell2mat(k);end;
if isfield(results,'flns_raw')
    results.b_raw(k) = results.flns_raw(k);
end
[~,k] = unique(results.ts,'first');
%results.b_flns = interp1(results.ts(k),results.flns(k),self.b_ts,'nearest');
%results.b_passindex= interp1(results.ts(k),results.passindex(k),self.b_ts,'nearest');

self.epoch = traj_epoch;

self.b_myvar = results.b_raw;
results.spk_raw = self.cel_myvar;

%self.b_myvar = results.b_flns;
results.spk_flns = self.cel_myvar;

%self.b_myvar = results.b_passindex;
results.spk_passindex = self.cel_myvar;

results.spk_theta = self.cel_theta;

if iscell(results.spk_passindex)
    results.spk_theta=cell2mat(results.spk_theta);
    results.spk_raw=cell2mat(results.spk_raw);
    results.spk_flns=cell2mat(results.spk_flns);
    results.spk_passindex=cell2mat(results.spk_passindex);
end
    results.spk_theta=mod(results.spk_theta,2*pi);
    
    [results.un, results.dn] = circcor(results.spk_theta,results.spk_passindex);

    
ts = self.cel_ts;
if ~iscell(ts)
    ts = {ts};
end
results.spk_theta_cycle = cellfun(@(x)interp1(results.lfpts,results.theta_cycle,x,'nearest'),ts,'UniformOutput',false);
if isscalar(results.spk_theta_cycle)
    results.spk_theta_cycle = results.spk_theta_cycle{1};
end


if ischar(siz)&&strcmpi(siz,'auto')
    n = floor(sqrt(length(plots)));
    m = n;
    nturn = false;
    while (n*m)<length(plots)
        if (nturn)
            n = n+1;
        else
            m = m+1;
        end
        nturn = ~nturn;
    end
else
    n=siz(1);
    m=siz(2);
end

if ischar(plot_inds)&&strcmpi(plot_inds,'auto')
    plot_inds = 1:length(plots);
end

for i=1:length(plots)
    subplot(n,m,plot_inds(i));
    if ~isempty(strfind(lower(plots{i}),'map'))&&(isempty(strfind(lower(plots{i}),':'))||all(strfind(lower(plots{i}),':')>strfind(lower(plots{i}),'map')))
        if ~isempty(strfind(lower(plots{i}),'polar'))||(isempty(strfind(lower(plots{i}),'spatial'))&&strcmpi(fieldness_func,'polar'))
            if ~isfield(results,'theta')
                varargin2 = varargin;
                if ~any(cellfun(@(x)isequal(x,'fieldness_func'),varargin2))
                    varargin2{end+1} = 'fieldness_func';
                end
                varargin2{find(cellfun(@(x)strcmpi(x,'fieldness_func'),varargin2))+1} = 'polar';
                resultshd = infieldness(self,cel,varargin2{:});
                
                results.theta = resultshd.theta;
                results.polar_rate_map = resultshd.rate_map;
                results.polar_flns_map = resultshd.flns_map;
            end
            
            if ~isfield(results,'polar_rate_map')
                if ~isempty(strfind(lower(plots{i}),'infieldness'))
                    rate_map = results.flns_map;
                elseif ~isempty(strfind(lower(plots{i}),'rate'))
                    rate_map = results.rate_map;
                end
            else
                if ~isempty(strfind(lower(plots{i}),'infieldness'))
                    rate_map = results.polar_flns_map;
                elseif ~isempty(strfind(lower(plots{i}),'rate'))
                    rate_map = results.polar_rate_map;
                end
            end
            
            theta = results.theta;
            
            
            h1=polar(deg2rad(theta(:)), rate_map(:), 'k'); hold on;
            set(h1,'linewidth',1.1);
            
            xs = rate_map(1:end-1).*cosd(theta(1:end-1)); % average
            ys = rate_map(1:end-1).*sind(theta(1:end-1));
            
            coordlims=axis;
            
            ang_hd = atan2(mean(ys),mean(xs)); % mean direction
            
            mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(rate_map(1:end-1)); % mean resultant length
            
            mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*coordlims(2); % for visualizations sake
            
            h1 = polar([ang_hd ang_hd], [0 mag_hd], 'r'); hold on
            set(h1,'linewidth',1.1)
            wu2 = self.HDWatsonsU2(cel);
            text(.55*xs(2), .8*ys(2), ['Watson''s U^2: ' num2str(wu2)], 'FontSize', 8);
            title(['Max, R:' num2str(mr,3)]);
            hold off;
        else
            if ~isfield(results,'xdim')
                varargin2 = varargin;
                if ~any(cellfun(@(x)isequal(x,'fieldness_func'),varargin2))
                    varargin2{end+1} = 'fieldness_func';
                end
                varargin2{find(cellfun(@(x)strcmpi(x,'fieldness_func'),varargin2))+1} = 'rate map';
                resultsspac = infieldness(self,cel,varargin2{:});
                
                results.xdim = resultsspac.xdim;
                results.ydim = resultsspac.ydim;
                results.spatial_rate_map = resultsspac.rate_map;
                results.spatial_flns_map = resultsspac.flns_map;
                results.no_occupancy = resultsspac.no_occupancy;
            end
            
            xdim = results.xdim;
            ydim = results.ydim;
            no_occupancy = results.no_occupancy;
            
            if ~isfield(results,'spatial_rate_map')
                if ~isempty(strfind(lower(plots{i}),'infieldness'))
                    unit = '';
                    rate_map = results.flns_map;
                elseif ~isempty(strfind(lower(plots{i}),'rate'))
                    unit = ' Hz';
                    rate_map = results.rate_map;
                end
            else
                if ~isempty(strfind(lower(plots{i}),'infieldness'))
                    unit = '';
                    rate_map = results.spatial_flns_map;
                elseif ~isempty(strfind(lower(plots{i}),'rate'))
                    unit = ' Hz';
                    [rate_map, xdim, ydim, ~, no_occupancy] = root.RateMap(root.cel);
                end
            end
            
            xdim = xdim(sum(~no_occupancy)~=0);
            ydim = ydim(sum(~no_occupancy')~=0);
            rate_map = rate_map(sum(~no_occupancy')~=0,sum(~no_occupancy)~=0);
            no_occupancy = no_occupancy(sum(~no_occupancy')~=0,sum(~no_occupancy)~=0);
            
            xs = [min(xdim) max(xdim)];
            ys = [min(ydim) max(ydim)];
            
            clims = [0 max(max(rate_map))];
            
            [cbar, clims] = CMBHOME.Utils.SmartColorbar(clims, 'jet(255)');
            
            rate_map(no_occupancy) = clims(1);
            
            imagesc(xdim, ydim, rate_map, clims); hold on;
            colormap(cbar);
            freezeColors;
            
            
            axis equal off
            
            xlim(diff(xs).*pad+xs);
            ylim(diff(ys).*pad+ys);
            
            line([xs(1)+.75*diff(xs), xs(2)], [-.03*diff(ys)+ys(1), -.03*diff(ys)+ys(1)], 'Color', 'k', 'LineWidth', 3);
            text(xs(2), -.03*diff(ys)+ys(1), [num2str(.25*diff(xs)*self.spatial_scale, 3) ' cm'], 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom');
            
            str_f = ['max: ' num2str(max(rate_map(:)), 2) unit ' mean: ' num2str(mean(rate_map(:)), 2) unit];
            
            text(xs(2), .045*diff(ys)+ys(2), str_f, 'FontSize', 6.8, 'FontWeight', 'bold', 'HorizontalAlign', 'right');
            
            set(gca,'YDir','normal'); % so plotting functions dont reverse axis
        end
        
    elseif ~isempty(strfind(lower(plots{i}),'trajectory'))
        
        hold on;
        
        self.epoch = infieldness_epoch;
        xs = self.x;
        ys = self.y;
        
        if iscell(self.x)
            cellfun(@(x,y)line(x,y,'Color',[1 1 1]*0.7,'LineWidth',1),self.x,self.y);
        else
            line(self.x,self.y,'Color',[1 1 1]*0.7,'LineWidth',1);
        end
        
        if ~isempty(strfind(lower(plots{i}),'rate map'))
            colors=jet(256);
            clims=[0 max(results.flns_raw)];
            spk_val = results.spk_raw;
        elseif ~isempty(strfind(lower(plots{i}),'pass index'))
            colors=jet(256);
            clims = [-1 1];
            spk_val = results.spk_passindex;
            colorbar('location','eastoutside');
        else
            colors=jet(255);
            spk_val = results.spk_flns;
            clims=[0 1];
        end
        
        if ~iscell(spk_val)
            spk_val = {spk_val};
        end
        
        self.epoch = traj_epoch;
        
        if iscell(self.x)
            cellfun(@(x,y)line(x,y,'Color','k','LineWidth',1),self.x,self.y);
            scatter(cell2mat(self.cel_x),cell2mat(self.cel_y), 10, cell2mat(spk_val), 'filled');
        else
            line(self.x,self.y,'Color','k','LineWidth',1);
            scatter(cell2mat(self.cel_x),cell2mat(self.cel_y), 10,cell2mat(spk_val),'filled');
        end
        
        colormap(colors);
        set(gca,'CLim',clims);
        
        hold off;
        
        axis equal
        
        axis off
        
        if ~iscell(xs), xs={xs};ys={ys}; end;
        xa = [min(cell2mat(xs)) max(cell2mat(xs))];
        ya= [min(cell2mat(ys)) max(cell2mat(ys))];
        
        xWindow = diff(xa).*pad+xa;
        yWindow = diff(ya).*pad+ya;
        
        xlim(xWindow);
        ylim(yWindow);
        
        line([xa(1)+.75*diff(xa), xa(2)], [-.03*diff(ya)+ya(1), -.03*diff(ya)+ya(1)], 'Color', 'k', 'LineWidth', 3);
        text(xa(2), -.03*diff(ya)+ya(1), [num2str(.25*diff(xa)*self.spatial_scale, 3) ' cm'], 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlign', 'right', 'VerticalAlign', 'bottom');
        
        %freezeColors;
    elseif ~isempty(strfind(lower(plots{i}),'scatter'))||~isempty(strfind(lower(plots{i}),'density'))
        
        self.epoch = traj_epoch;
        p = inputParser;
        p.addOptional('x',results.spk_passindex);
        p.addOptional('y',mod(results.spk_theta,2*pi));
        p.addParamValue('repeat',[0 0;0 1]);
        p.addParamValue('range',[-1 1;0 360]);
        p.addParamValue('xLabel','Pass index');
        p.addParamValue('yLabel','Theta phase (deg)');
        k = eval(['{' plots{i}(strfind(plots{i},':')+1:end) '}']);
        p.parse(k{:});
        
        x = p.Results.x;
        y = p.Results.y;
        
        repeat = p.Results.repeat;
        rng = p.Results.range;
        xLabel = p.Results.xLabel;
        yLabel = p.Results.yLabel;
        
        if iscell(x), x=cell2mat(x); end;
        if iscell(y), y=cell2mat(y); end;
        
        if isempty(strfind(lower(plots{i}),':'))
            y = rad2deg(y);
        end
        
        if any(repeat(1,:))
            x = repmat(x,[range(repeat(1,:))+1 1])+reshape(repmat(repeat(1,1):repeat(1,2),[1 numel(x)]),[numel(x)*(range(repeat(1,:))+1) 1])*range(rng(1,:));
            y = repmat(y,[range(repeat(1,:))+1 1]);
        end
        
        if any(repeat(2,:))
            y = repmat(y,[range(repeat(2,:))+1 1])+reshape(repmat((repeat(2,1):repeat(2,2))',[1 numel(y)])',[numel(y)*(range(repeat(2,:))+1) 1])*range(rng(2,:));
            x = repmat(x,[range(repeat(2,:))+1 1]);
        end
        
        if ~isempty(strfind(lower(plots{i}),'scatter'))
            scatter(x,y,5,'o','filled');
            
            xlim(range(rng(1,:))*repeat(1,:)+rng(1,:));
            ylim(range(rng(2,:))*repeat(2,:)+rng(2,:));
            
            k = ylim;
            [rho,p,s]=anglecor2(results.spk_passindex,results.spk_theta);
results.rho = rho;
results.p = p;
results.s = s;
results.is_precessing = p<0.05&&rad2deg(s*2)<-22&&rad2deg(s*2)>-1440;
%             keyboard
            title(['rho=' num2str(rho,'%11.2g') 'p=' num2str(p,'%11.2g') 's=' num2str(rad2deg(s*2),'%11.2g') 'deg/pass']);
        else
            
%             keyboard
            %%
%            self.b_myvar=results.b_passindex;
            ns = histcn([rad2deg(results.spk_theta) results.spk_passindex],linspace(0,360,101),linspace(-1,1,41));
            
            t = cell2mat2(self.ts);
            pii = cell2mat2(self.myvar);
            theta = cell2mat2(self.lfp.theta_phase);
            ts = cell2mat2(self.lfp.ts);
            
            [~,temp]=unique(t);
            [~,temp2]=unique(ts);
%             
%             keyboard
            oc = histcn([rad2deg(mod(theta(temp2),2*pi)) nanInterp1(t(temp),pii(temp),unique(ts),'nearest')],linspace(0,360,101),linspace(-1,1,41))*mean(diff(ts));
            
            rm = ns./oc;
            rm(isnan(rm))=0;
            
            h = fspecial('gaussian',3,1);
            rm = imfilter(rm,h,'replicate');
            
            imagesc(linspace(0,1,40),linspace(0,360*2,100),repmat(rm,[2 1]));
            set(gca,'YDir','Normal');
            xlabel('Pass index');ylabel('Theta phase (deg)');
            title(['max: ' num2str(max(rm(:)),'%2.2f') ' Hz']);
            freezeColors;
            %%
        end
        
        
        
        xlabel(xLabel);ylabel(yLabel);
        axis square;
    elseif strcmpi(plots{i},'Histogram:Pass Index')
        if iscell(results.spk_passindex)
            [nspike,xaxis] = hist(cell2mat(results.spk_passindex),20);
        else [nspike,xaxis] = hist(results.spk_passindex,20); end;
        
        
%        self.b_myvar = results.b_passindex;
        self.epoch = traj_epoch;
        
        if iscell(self.myvar)
            bintime=hist(cell2mat(self.myvar),20);
        else
            bintime=hist(self.myvar,20);
        end
        
        bar(xaxis,nspike./bintime*self.fs_video,1);
        xlim([-1 1]);
        xlabel('Field Pass Index');ylabel('Spike rate (Hz)');
    elseif strcmpi(plots{i},'Histogram:Theta')
%         keyboard
        %%
        self.epoch = traj_epoch;
        
        if iscell(results.spk_theta)
            [nspike,xaxis] = hist(rad2deg(cell2mat(results.spk_theta)),20);
        else [nspike,xaxis] = hist(rad2deg(results.spk_theta),20); end;
        
        if iscell(self.lfp.theta_phase)
            bintime=hist(rad2deg(cell2mat(self.lfp.theta_phase)),20);
        else
            bintime=hist(rad2deg(self.lfp.theta_phase),20);
        end
        
        bar([xaxis xaxis+360],repmat(nspike./bintime*self.lfp.fs,[1 2]),1);
        xlim([0 360*2]);
        xlabel('Theta Phase (deg)');ylabel('Spike rate (Hz)');
        
    elseif strcmpi(plots{i},'Correl')
        %%
        self.epoch = traj_epoch;
        
        [sf,theta_index]=self.IntrinsicFrequency2(self.cel);
        thet = self.lfp.theta;
        if iscell(thet), thet=cell2mat(thet); end;
        
        pxx=periodogram(thet);
        
        f = length(pxx);
        inds=(round(6/(self.lfp.fs/2)*length(pxx)):round(10/(self.lfp.fs/2)*length(pxx)))';
        pxx = pxx(inds);
        f=sum(inds.*pxx)/sum(pxx)/f*(self.lfp.fs/2);
        
        k = xlim;
        xlim([0 k(2)]);
        xaxis = 0:0.01:k(2);
        k = ylim;
        ylim([0 k(2)]);
        plot(xaxis,k(2)/4*cos(f*(2*pi)*xaxis)+k(2)/2,'b--','LineWidth',2);
        legend off;
        text('string',['{\color{Red} Fit:' num2str(sf,'%11.3g') 'Hz}    {\color{Blue} EEG Average Frequency:' num2str(f,'%11.3g') 'HZ}\newline{\color{Red}Theta index:' num2str(theta_index,'%11.3g') '}'],'position',[0.05 0.9*k(2)],'interpreter','tex')
        
    end
    
    
    %pause(1/100);
end

end
