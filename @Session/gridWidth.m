function gw = gridWidth(self,cel)
% gridWidth = root.gridWidth(root.cel)
%
% Returns the diameter of the grid field, based on the center peak in the
% spatial autocorrelogram. Width is in PIXELS (divide by 2 to get spatial
% rather than autocorr spatial)
%
% Methods from:
% Anisotropic encoding of three-dimensional space by place cells and grid
% cells. Heyman et al 2012. Nature Neuroscience.  

% wchapman 20140404

    rate_map = self.RateMap(cel, 'continuize_epochs', 1);
    ac = CMBHOME.Utils.moserac(rate_map);
    
    ac = ac>0.15;
    
    rp = regionprops(ac,'MajorAxisLength','centroid','PixelIdxList');

    mp = floor(numel(ac)/2);    % middle of the ac
    for i = 1:size(rp,1)
        if ismember(mp,rp(i).PixelIdxList)
            ind = i;
        end
    end
    
    gw = rp(ind).MajorAxisLength;
    
end