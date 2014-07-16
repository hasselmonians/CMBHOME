function [S, f] = plot_spk_spectrum(self, cel)
% root.plot_spk_spectrum(cel);
%
% plots the Frequency power spectrum for spiking of cell cel for every
% root.epoch.
%
% Returns
% S -> spectrum
% F -> frequency vector
%
% andrew bogaard 3 april 2010

    import CMBHOME.Utils.*
    
    AddChronuxPackage;

    params.Fs = 40; % set params for spect function. bin in .025s bins
    params.fpass = [0 20];
    params.pad = 6;
    
    spk_ts = self.spk_ts(cel);
    
    if iscell(spk_ts)
        if length(spk_ts)>10
            reply = input('You requested to plot more than 10 epochs. This could take awhile. Continue? Y/N [Y]: ', 's');
            if isempty(reply) || strcmpi(reply, 'n')
                return
            end
            
            S = cell(length(spk_ts,1);
            f = cell(length(spk_ts,1);
    
        end
        
        for i =1:length(spk_ts)
            
            ts = histc(spk_ts{i}, self.epoch(i,1):params.Fs^-1:self.epoch(i,2));
            
            params.tapers = [ .0013 diff(self.epoch(i,:)) 1 ];
            
            [S{i}, f{i}] = mtspectrumpb(ts, params); % there exists a point time folder (so that we dont have to bin), but I find that it doesnt work well
            
            %line(f{i}, 10*log10(S{i}), 'Color', 'k', 'LineWidth', 1), hold on;
            line(f{i}, S{i}, 'Color', 'k', 'LineWidth', 1), hold on;
                        
        end
    else
        
        ts = histc(spk_ts, self.epoch(1):params.Fs^-1:self.epoch(2));
        
        params.tapers = [ .0013 diff(self.epoch) 1 ];
        
        [S, f] = mtspectrumpb(ts, params);
        
        %line(f, 10*log10(S), 'Color', 'k', 'LineWidth', 1);
        line(f, S, 'Color', 'k', 'LineWidth', 1);
        
    end

    
end
