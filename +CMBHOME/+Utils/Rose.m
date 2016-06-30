function [mra mrl r t ul ll] = Rose(xin,be,ifPlot,ifSmooth)
% [mra mrl r t ul ll] = rose2(alpha,binEdges,ifPlot,ifSmooth)
    if ~exist('be','var')
        be = [];
    end
    
    if isempty(be)
        be = linspace(0,2*pi,16);
    end
    
    if ~exist('ifSmooth','var')
        ifSmooth = 0;
    end

    xin(isnan(xin)) = [];
    xin(isinf(xin)) = [];
    
    %# data and angular histogram
    [t,r] = rose(xin,be);
    r = r(:); t =t(:);
    mrl = CMBHOME.Utils.circ.circ_r(xin);
    mra = CMBHOME.Utils.circ.circ_mean(xin);
 
    
    %% threshes
    t1 = t(2:4:end);
    t2 = t(3:4:end);
    tm = mean([t1(:) t2(:)],2);
    
    r1 = r(2:4:end);
    r2 = r(3:4:end);
    rm = mean([r2(:) r2(:)],2);
    
    badt = tm(find(rm<(max(rm)/3)));
    mra = CMBHOME.Utils.circ.circ_mean(xin);
    dv = min(wrapTo2Pi(mra - badt));
    
    ul = mra+dv;
    ll = mra-dv;
    
    
    %% Plot
    if ~exist('ifPlot','var')
        ifPlot = 1;
    end
    
    if ifPlot==1
        %# set plot's max radial ticks
        %figure
        rMax = max(r);
        %rMax = 50*(ceil(rMax/50));
        h = polar(0, rMax);
        delete(h)
        set(gca, 'Nextplot','add')
        
        %# draw patches instead of lines: polar(t,r)
        [x,y] = pol2cart(t,r);
        h = patch(reshape(x,4,[]), reshape(y,4,[]), [0.7 0.7 0.7]);
        %alpha(h, 0.5)       %# note: this switches to OpenGL renderer
        set(gcf,'Renderer','painters')
        
        % put the mean-resultant on the plot as well
        [x y] = pol2cart(mra,mrl*rMax);
        plot([0 x],[0 y],'r','LineWidth',2)
    end
    
end
