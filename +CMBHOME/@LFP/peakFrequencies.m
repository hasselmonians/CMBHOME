function [peakTheta, peakLowGamma, peakHighGamma] = peakFrequencies(lfp)
%Determines the peak frequencies in theta, low gammma, and high gamma bands
%of significant theta epochs in an lfp.

%check to see if valid LFP channel
if ~all(lfp.signal == 0)
    %determine theta epochs
    [thetaEpochs, ~, ~] = ThetaEpochs(lfp);

    if ~isempty(thetaEpochs)
        %determine power spectrum for each epoch
        [r, ~] = size(thetaEpochs);
        out = struct();

        for i = 1:r
            signal = lfp.signal(find(lfp.ts == thetaEpochs(i,1)) : find(lfp.ts == thetaEpochs(i,2)));
            L = length(signal);
            NFFT = 2^nextpow2(L);
            Y = fft(signal, NFFT) / L;
            out(i).f = lfp.fs / 2*linspace(0, 1, NFFT/2+1);
            out(i).FFT = 2*abs(Y(1:NFFT/2+1));
        end

        %find minimum FFT length
        lens = nan(1,length(out));
        for i = 1:length(out)
            lens(i) = length(out(i).f);
        end
        minLen = min(lens);

        %downsample longer FFTs
        for i = 1:length(out)
            n = round(length(out(i).f) / minLen);
            if n > 1
                out(i).f = out(i).f(1:n:n*minLen-(n-1));
                out(i).FFT = out(i).FFT(1:n:n*minLen-(n-1));
            end
        end

        %average power spectra
        FFTs = zeros(r,minLen);
        for i = 1:r
            FFTs(i,:) = out(i).FFT;
        end
        f = out(1).f;
        FFTavg = mean(FFTs);

        %determine peak frequencies
        theta = find(f >= 4 & f <= 12);
        lowGamma = find(f >= 20 & f <= 60);
        highGamma = find(f >= 65 & f <= 150);

        thetaInd = theta(1) + find(FFTavg(theta) == max(FFTavg(theta))) - 1;
        lowGammaInd = lowGamma(1) + find(FFTavg(lowGamma) == max(FFTavg(lowGamma))) - 1;
        highGammaInd = highGamma(1) + find(FFTavg(highGamma) == max(FFTavg(highGamma))) - 1;

        peakTheta = f(thetaInd);
        peakLowGamma = f(lowGammaInd);
        peakHighGamma = f(highGammaInd);
    else
        peakTheta = NaN;
        peakLowGamma = NaN;
        peakHighGamma = NaN;
    end
    
else
    peakTheta = NaN;
    peakLowGamma = NaN;
    peakHighGamma = NaN;
end

end