function root = AlignSpike2Session(root)
% root = root.AppendThetatoSpike;
%
% For every spike object in root.spike, realigns spike times to root.b_ts,
% and deletes the rest of the cells information other than label
%
% andrew 23 sept 2011

for j = 1:numel(root.spike)
    if ~isempty(root.spike(j).ts)
        root.spike(j) = CMBHOME.Spike('ts', root.spike(j).ts,...
                                'vid_ts', root.b_ts,...
                                'tet_name', root.spike(j).tet_name ,...
                                'prop', root.spike(j).prop);
    end

end
    