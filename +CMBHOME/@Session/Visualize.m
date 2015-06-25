function Visualize(self, varargin)

% Prints an 8.5x11 figure of subplots indicated by varargin. See options
% below. Optionally prints to pdf

import CMBHOME.Utils.*

p = inputParser;

p.addRequired('self', @(x) all(class(x)=='CMBHOME.Session')); % must be CMB Standard Session class

p.addParamValue('cells',            [], @(x) size(x,2)==2); % Nx2 matrix of N cells like [tetrode, cell;...]
p.addParamValue('screen_print',     'off' , @(x) strcmp(x, 'off') | strcmp(x,'on')); % boolean indicating whether or not to print figure to screen and pause at each
p.addParamValue('header',           1 , @(x) length(x)==1); % done
p.addParamValue('rate_map',         0 , @(x) length(x)==1); % done
p.addParamValue('polar_rate_map',   0 , @(x) length(x)==1); % done
p.addParamValue('auto_corr',        0 , @(x) length(x)==1); % done
p.addParamValue('head_dir',         0 , @(x) length(x)==1); % done
p.addParamValue('trajectory',       0 , @(x) length(x)==1); % done
p.addParamValue('isi_dist',         0 , @(x) length(x)==1); % done
p.addParamValue('theta_phase_dist', 0 , @(x) length(x)==1); % done
p.addParamValue('ff_vel',           0 , @(x) length(x)==1);
p.addParamValue('raster',           [], @(x) size(x, 2)==3); % if raster plot, submit [start time, stop time, offset]
p.addParamValue('pdf_fname',        [], @(x) ischar(x)); % string name of pdf file, like 'myPDF.pdf'

p.parse(self, varargin{:});

PDF_fname = p.Results.pdf_fname;
cells = p.Results.cells;
screen_print = p.Results.screen_print;

fig = figure('Position', [50, -150, 650, 1100], 'PaperPosition', [1 .2 6.5 10.6], 'Visible', screen_print, 'Color', 'w'); 

if isempty(self.spike), disp('No Cells'); return; end

if isempty(cells), cells = self.cells; end

for j = 1:size(self.cells,1)
   
    fig_ind = 1; % start at the top of a new subplot
    
    if p.Results.header
        figure(fig), subplot(4, 2, fig_ind),
        f_header(self, cells(j,:));
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if p.Results.trajectory
        figure(fig), subplot(4, 2, fig_ind),
        f_trajectory(self, cells(j,:));
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
        
    if p.Results.rate_map
        figure(fig), subplot(4, 2, fig_ind),
        f_rate_map(self, cells(j,:));
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if p.Results.polar_rate_map
        figure(fig), subplot(4, 2, fig_ind),
        f_polar_rate_map(self, cells(j,:));
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if p.Results.auto_corr
        figure(fig), subplot(4, 2, fig_ind),
        f_auto_corr(self, cells(j,:), 2, .150);
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if p.Results.head_dir
        figure(fig), subplot(4, 2, fig_ind),
        f_head_dir(self);
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if p.Results.isi_dist
        figure(fig), subplot(4, 2, fig_ind),
        f_isi_dist(self, cells(j,:));
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if p.Results.theta_phase_dist
        figure(fig), subplot(4, 2, fig_ind),
        f_theta_phase_dist(self, cells(j,:));
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if p.Results.ff_vel
        figure(fig), subplot(4, 2, fig_ind),
        f_ff_vel(self, cells(j,:));
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
        
    if p.Results.raster
        figure(fig), subplot(4, 2, fig_ind),
        f_raster(self, cells(j,:), raster);
        fig_ind = IncrementSP(fig_ind, PDF_fname);
    end
    
    if strcmp(screen_print, 'on'), pause; end
           
    if PDF_fname, print('-dpsc2',[PDF_fname '.ps'], '-append', ['-f' num2str(gcf)]); end % append pdf ps file

    clf;
   
end

if PDF_fname, ps2pdf('psfile',[PDF_fname '.ps'], 'pdffile', PDF_fname, 'deletepsfile',1); end % convert the ps to pdf 
end

function i = IncrementSP(i, PDF_fname)
i = i+1;
       
if i > 8

    i = 1;

    if PDF_fname, print('-dpsc2',[PDF_fname '.ps'], '-append', ['-f' num2str(gcf)]); end % append pdf ps file

    clf;

end
end


