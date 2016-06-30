function root = alpha2counting(root)

%Takes in an CMBHOME object and reorders the root.pathlfp to be in counting
%order (eg: 'DATE_BASE.EEG1', 'DATE_BASE.EEG2')

%Bill 2012.05.23

for i = 1:length(root.path_lfp)
    eegs.names(i) = root.path_lfp{i}(3);
    eegs.lengths(i)= length(eegs.names{i});
end

[trash ind] =sort(eegs.lengths);

for i = 1:length(ind)
    root.path_lfp{i}(3)=eegs.names(ind(i));
end



