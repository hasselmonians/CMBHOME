function plot_lfp_spectrum(self, TBP, n_tapers, f_range, theta_analysis)
% root.plot_lfp_spectrogram(lfp_ind, windowsize, windowsin, bandwidth,
% f_range);
%
% Plots a spectrogram in jet colormap for root.b_lfp(lfp_ind). Params for
% the moving FFT analysis above.
%
% Uses the Chronux toolbox.
%
% andrew bogaard 3 april 2010

    import CMBHOME.Utils.*
    
    %AddChronuxPackage;
       
    if ~exist('TBP', 'var')
        TBP = 10;
    end
    
    if ~exist('n_tapers', 'var')
        n_tapers = 9;
    end
    
    if ~exist('theta_analysis', 'var')
        theta_analysis = 0;
    end
        
    signal = self.lfp.signal;
    
    if iscell(signal), signal = vertcat(signal{:}); end
        
    params.tapers = [TBP, n_tapers];
    params.fpass = f_range;
    params.Fs = self.lfp.fs;
    params.pad = 2;

    [S,f]=mtspectrumc(signal,params);

    stdf = .2/mean(diff(f)); % elements in a .2 hz std
    
    kernel = pdf('Normal', -3*stdf:3*stdf, 0, stdf); % smooth with gaussian kernel
    
    S = conv(S, kernel, 'same'); % plot(f, S)
    
    [xmax,imax] = extrema(S); % find local maxima and see if any lie within theta

    ftmp = f(imax);
    stmp = S(imax);

    [a_peak,ind] = max(stmp(ftmp > 7 & ftmp < 11)); % intrinsic frequency

    ftmp = ftmp(ftmp > 7 & ftmp < 11);

    if ~isempty(ind)                                    % if theta peak was found
        f_peak = ftmp(ind)';     
        a_peak_av = mean(S(f>f_peak-1 & f<f_peak+1));
        peak = S(f==f_peak);
    else                                                % if theta peak was not found
        f_peak = 0;
        peak = 0;
        a_peak_av = 0;
    end

    power_ratio = a_peak_av / mean(S);
        
    %line(f, 10*log10(S)', 'Color', 'k', 'LineWidth', 1.5)
    line(f, S', 'Color', 'k', 'LineWidth', 1.5), hold on
    
    a_peak_av = mean(S(f>f_peak-1 & f<f_peak+1));

    signal_power = mean(signal.^2);

    if a_peak_av<1.5*signal_power
        f_peak = [];
    end
    
    ys = ylim;
    if theta_analysis
        line([f_peak f_peak], ys, 'Color', 'r', 'linewidth', 1.5), hold on;
        text(f_peak+15*.02, ys(1)+diff(ys)*.02, ['Intrinsic Frequency = ' num2str(f_peak) ' Hz'], 'FontSize', 8);
        text(f_peak+15*.02, ys(1)+diff(ys)*.015, ['Power_theta / Power_av = ' num2str(power_ratio)], 'FontSize', 8);
    end

end
