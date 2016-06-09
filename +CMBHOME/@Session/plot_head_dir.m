function plot_head_dir(self)
    % Plots the head direction vector on the current axis
    line(CMBHOME.Utils.ContinuizeEpochs(self.ts), CMBHOME.Utils.ContinuizeEpochs(self.headdir), 'Color', 'k');
end

