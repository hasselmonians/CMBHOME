function root = AddLFPTheta(root)
% Adds filtered LFP signal (and phase, amplitude) to root.b_lfp
%
% For all root.b_lfp, adds theta filtered signal and theta phase vectors to
% speed up further analysis

if isempty(root.b_lfp)
    
    disp('No LFP data was loaded, adding all root.path_lfp');
    
    root = root.LoadLFP(1:length(root.path_lfp));
    
    if isempty(root.b_lfp), error('No LFP files were loaded'); end
    
end

for i = 1:length(root.b_lfp)

    root.b_lfp(i) = root.b_lfp(i).AppendTheta;
    
    root.b_lfp(i) = root.b_lfp(i).AppendThetaPhase;
    
    root.b_lfp(i) = root.b_lfp(i).AppendThetaAmplitude;
    
end