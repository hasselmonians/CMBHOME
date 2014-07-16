function [tuning_curve, theta, ang_hd, mr] = DirectionalTuningFcn(self, cel, varargin)
% (1) [tuning_curve, theta] = root.DirectionalTuningFcn(cel);
% (2) [tuning_curve, theta] = root.DirectionalTuningFcn(cel, 'binsize', degrees, 'Continuize', 1);
% 
% Computes the directional (orientational) tuning curve. Also called rate
% map, 
% 
% Arguments (see syntax (2)):
%
% binsize -> width of binsize in degrees
% Continuize -> merge overlapping epochs, and concatonate all data from
% root.epoch
%
% Returns
%
% tuning_curve -> NxMxO matrix where N is the number of bins as per
% 'binsize', O is the number of epochs and M is the number of cells
%
% theta -> vector of bin centers in degrees


    import CMBHOME.Utils.*
    
    p = inputParser;

    p.addRequired('self');
    p.addRequired('cel');
    p.addParamValue('binsize', 6, @(x) isnumeric(x)); % 6 degrees default binsize
    p.addParamValue('Continuize',  0, @(x) numel(x)==1);

    p.parse(self, cel, varargin{:});

    binsize = p.Results.binsize;
    Continuize = p.Results.Continuize;

    theta=-180+binsize/2:binsize:180-binsize/2;
    
    theta = theta(:);

    angle_occupancy = self.DirectionalOccupancy(binsize, Continuize); % bin x 1 epoch cell array of vectors of seconds at theta. if continuize, epoch = 1
        
    if Continuize && size(self.epoch,1)>1 % continuize multiple epochs
        
        self = MergeEpochs(self);
        
        spk_headdir = ContinuizeEpochs(self.spk_headdir(cel));

    elseif size(cel,1)==1 && size(self.epoch,1)==1 % if one cell, one epoch
        
        spk_headdir = self.spk_headdir(cel);      
        
    elseif size(self.epoch,1)>1 || size(cel,1)>1 % if multiple epochs and not continuized, or multiple cells
        
        spk_headdir = self.spk_headdir(cel); % cell array of epochs x cells

        spk_headdir = CatCA(spk_headdir); % function in CMBHOME.Utils converts cell array to matrix

    end
    
    num_spikes = zeros(length(theta), size(cel, 1), size(spk_headdir,3));
    
    for i = 1:size(spk_headdir, 3) % if multiple cells and epochs. this partially vectorizes the functions since hist accepts matrices
        
        num_spikes(:,:,i) = hist(spk_headdir(:,:,i), theta); 
        
    end

    if Continuize
        tuning_curve = num_spikes ./ repmat(angle_occupancy, [1 size(cel,1), size(spk_headdir,3)]);
    else
        tuning_curve = num_spikes ./ repmat(angle_occupancy, [1 size(cel,1)]);
    end
    
    [ang_hd, mr] = GetOtherStats(tuning_curve, theta); 
    
    
end

function [ang_hd, mr] = GetOtherStats(tuning_curve, theta)

    tf = ~isnan(tuning_curve); % where we arent nan

    theta=theta*unitsratio('rad','deg');
    
    theta = theta(tf); % remove nans
    
    tuning_curve = tuning_curve(tf);

    xs = tuning_curve.*cos(theta); % average 
    ys = tuning_curve.*sin(theta);

    ang_hd = atan2(mean(ys),mean(xs)); % mean direction
    
    mr = (cos(ang_hd)*sum(xs) + sin(ang_hd)*sum(ys)) / sum(tuning_curve); % mean resultant length
    
end