function auto_corr_rm = plot_rate_map_ac(self, cel, rate_map, xs, ys)
% Peristimulus time histogram of spikes between times

    import CMBHOME.Utils.*
        
    if ~exist('rate_map', 'var'), rate_map=[]; end
    
    if ~isempty(rate_map) %&& exist('sx', 'var') && exist('ys', 'var')
               
        %rate_map = rate_map - mean(mean(rate_map));
        
        auto_corr_rm = moserac(rate_map);
        
%         auto_corr_rm = xcorr2(rate_map);
%         
%         [ c, r ] = meshgrid(1:size(auto_corr_rm,2), 1:size(auto_corr_rm,1));
%         
%         c = min(c, abs(c-size(auto_corr_rm,2)-1));
%         
%         r = min(r, abs(r-size(auto_corr_rm,1)-1));
%         
%         norm = c.*r;
%         
%         cc = gcf; 
%         
%         figure, imagesc(norm), 
%         
%         figure(cc)
%       
%         auto_corr_rm = auto_corr_rm./norm;
        
    else
        
        [rate_map, xs, ys] = self.RateMap(cel, 'continuize_epochs', 1);
        
        %rate_map = rate_map - mean(mean(rate_map));
        
        auto_corr_rm = moserac(rate_map);
        
    end
       
    clims = [min(auto_corr_rm(:)), .8];
    
    [cbar, clims] = CMBHOME.Utils.SmartColorbar(clims, 'jet(255)');

    auto_corr_rm(isnan(auto_corr_rm)) = clims(1);
    
    if any(isnan(clims)), text(.01, .3, 'No figure', 'Fontsize', 8), axis off, return; 
    else if ~diff(clims)>0, text(.01, .3, 'No figure', 'Fontsize', 8), axis off, return; end
    
    imagesc(auto_corr_rm, clims);
    
    colormap(cbar);
    
    %set(t, 'AlphaDataMapping', 'none');   %   until this works properly
                                          %without screwing up other axes in the subplot, screw it!
    %set(gca, 'DrawMode', 'fast');
    %set(t,'AlphaData', ~isnan(auto_corr_rm));

    axis equal
    
    axis off

    set(gca, 'Box', 'on')
    
    title(['Rate Map ACorr T' int2str(cel(1)) 'C' int2str(cel(2))]);
    
   % ylabel('pixels');
   % xlabel('pixels');
    
    [ymax, xmax] = size(auto_corr_rm); % set axis limits
    
    xlim([1, xmax]);
    ylim([1, ymax]);
    
end
