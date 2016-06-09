function root = AlignSpike2Session(root)
% Aligns Spike times to the Session frames. Populates root.Spike(n,m).i

%
% For every spike object in root.spike, realigns spike times to root.b_ts,
% and deletes the rest of the cells information other than label

for j = 1:numel(root.spike)
    if ~isempty(root.spike(j).ts)
        root.spike(j) = CMBHOME.Spike('ts', root.spike(j).ts,...
                                'vid_ts', root.b_ts,...
                                'tet_name', root.spike(j).tet_name ,...
                                'prop', root.spike(j).prop);
    end

end
    