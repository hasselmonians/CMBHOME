function plot_isi_dist(self, cel)

    ts = self.spk_ts(cel);
    
    if iscell(ts)
        
        tmp_ts = cell(length(unique(self.b_epoch_group)), 1); % make cell array of ISIs

        ind = 1;

        for i = unique(self.b_epoch_group)'

            dts = cellfun(@diff, ts(self.b_epoch_group==i), 'Uni', false);

            tmp_ts{ind} = vertcat(dts{:});

            ind = ind+1;

        end

        ts = tmp_ts;
        
        if length(ts)<4
        
            colors = cool(length(ts)).*.75;

            for i = 1:length(ts)

                 hist(log10(ts{i}),100); 

                 hold on;

                 h = findobj(gca,'Type','patch');

                 set(h(1),'FaceColor',colors(i,:),'EdgeColor',colors(i,:), 'FaceAlpha', .4, 'EdgeAlpha', .4);
             
            end
            
            h_legend=legend(self.active_event{:,1});
            set(h_legend,'FontSize',7);
            
        else
            
            ts = vertcat(ts{:});
            
            hist(log10(ts),100);
        
            h = findobj(gca,'Type','patch');

            set(h,'FaceColor','k','EdgeColor','k')
            
            h_legend=legend('Combined Epochs');
            set(h_legend,'FontSize',7);
            
        end
        
    else
        
        d_ts = diff(ts);
        
        hist(log10(d_ts),100);
        
        h = findobj(gca,'Type','patch');
    
        set(h,'FaceColor','k','EdgeColor','k')
        
        h_legend=legend(self.active_event{:,1});
        set(h_legend,'FontSize',7);
        
    end
    
    title('ISI Distribution');

    set(gca,'xlim',[-4 4],'xtick',-3:3,'xticklabel',{'0.001','0.01','0.1','1','10','100','1000'});

    xlabel('Interspike Interval (seconds)'); ylabel('Count');
    

        
end

