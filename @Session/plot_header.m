function plot_header(self, cel)
% root.plot_header(cel)
%
% prints the tetrode and cell number, current event labels, average firing
% frequency for these epochs.
%
% Used in Visualize2 primarily.
%
% andrew 3 april 2010

    % find event strings
    if all(arrayfun(@(a) any(a==self.event{:,2}), self.epoch)) % if both timestamps exist in structin.events
        eventInds = arrayfun(@(a) find(a==self.event), self.epoch);
        str1 = self.event{eventInds(1),1};
        str2 = self.event{eventInds(2),1};
    else
        str1 = ['Unknown event label. t_s_t_a_r_t = ' num2str(self.epoch(1))];
        str2 = ['Unknown event label. t_s_t_o_p = ' num2str(self.epoch(2))];
    end
    axis off
    text(0, .95, ['Tetrode ' num2str(cel(1))], 'FontSize', 14);
    text(0, .7, ['Cell ' num2str(cel(2))], 'FontSize', 14);
    text(0, .4, ['Event Start: ' str1], 'FontSize', 8);
    text(0, .25, ['Event Stop: ' str2], 'FontSize', 8);

    N_spikes = length(self.spk_ts(cel(1), cel(2)));

    F = N_spikes / diff(self.epoch);

    text(0, .1, ['Spikes/Sec during this event ' num2str(F)], 'FontSize', 8);

end

