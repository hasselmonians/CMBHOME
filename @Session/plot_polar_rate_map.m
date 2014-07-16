function [ratemap, ang_hd, mr, peak] = plot_polar_rate_map(self, cel, params, MaxRho)
% root.plot_polar_rate_map(cel);
% root.plot_polar_rate_map(cel, params, MaxRho);
%
% params is optional - [title, u2 score] (default [1 1])
% MaxRho is optional - positive number, upper bound of Rho
%
% andrew 3 april 2010

if ~exist('params', 'var'), params = [1 1]; end

    import CMBHOME.Utils.*

    binsize=6; %binsize in degrees
       
    [ratemap, theta] = self.DirectionalTuningFcn(cel, 'binsize', binsize, 'Continuize', 1);
    
    theta=theta*unitsratio('rad','deg');
    
    theta = cat(1, theta(:), theta(1));
    
    ratemap = cat(1, ratemap(:), ratemap(1));
    
    peak = max(ratemap);
    
    if ~exist('MaxRho', 'var')
        h1=polar(theta(:), ratemap(:), 'k'); hold on;
    else
        if MaxRho=='yes'
            h1=myPolar(theta(:), ratemap(:), max(ratemap), 'k'); hold on;
        else
            h1=myPolar(theta(:), ratemap(:), MaxRho, 'k'); hold on;
        end
    end
    
    set(h1,'linewidth',1.1)

    xs = ratemap(1:end-1).*cos(theta(1:end-1)); % average 
    ys = ratemap(1:end-1).*sin(theta(1:end-1));

    coordlims=axis;

    ang_hd = atan2(mean(ys),mean(xs)); % mean direction
    
    mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(ratemap(1:end-1)); % mean resultant length
    
    mag_hd = sqrt(sum(ys)^2+sum(xs)^2)/sqrt(sum(abs(ys))^2+sum(abs(xs))^2)*coordlims(2); % for visualizations sake

    if ~exist('MaxRho', 'var')
        h1 = polar([ang_hd ang_hd], [0 mag_hd], 'r'); hold on
    else
        if MaxRho=='yes'
            h1=myPolar(theta(:), ratemap(:), max(ratemap), 'k'); hold on;
        else
            h1=myPolar(theta(:), ratemap(:), MaxRho, 'k'); hold on;
        end
    end
    
    set(h1,'linewidth',1.1)

    xs = xlim;
    ys = ylim;

    if params(2)    
        wu2 = self.HDWatsonsU2(cel); % calculate watsons u2
        if numel(wu2)==1  
            text(.55*xs(2), .8*ys(2), ['Watson''s U^2: ' num2str(wu2)], 'FontSize', 8); 
        end
    end
    
    if params(1), title(['Head Direction Firing Rate, R:' num2str(mr,3)]); end
    
end