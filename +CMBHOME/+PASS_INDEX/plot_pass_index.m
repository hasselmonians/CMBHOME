function [ results ] = plot_pass_index( varargin )

import CMBHOME.PASS_INDEX.*;

ip = plot_pass_index_parser(varargin{:});
p = pass_index_parser(varargin{:});

for i = fields(p.Results)'
    eval([i{1} ' = p.Results.' i{1} ';']);
end
for i = fields(ip.Results)'
    eval([i{1} ' = ip.Results.' i{1} ';']);
end
for i = fields(results)'
    eval([i{1} ' = results.' i{1} ';']);
end

clear ip;

set(gcf,'color','w');
for i=1:numel(plots)
    subplot(subplots{i,:});
    switch plots{i}
        case 'Trajectory'
            switch size(pos,2)
                case 1 % 1D
                    plot(pos_ts,pos);
                    xlim(minmax(pos_ts(:)'));
                    ylim([-1 1]*range(pos(:))*1.1+nanmean(pos(:)));
                    
                    xlabel('Time (s)');
                    ylabel(['State (' units ')']);
                    
                    spkpos = spk_pos(pos_ts,pos,spk_ts);
                    hold on;
                    scatter(spk_ts,spkpos,'o','filled','CData',spk_pass_index,'SizeData',20);
                    hold off
                case 2 %2D
                    plot(pos(:,1),pos(:,2),'k');
                    spkpos = spk_pos(pos_ts,pos,spk_ts);
                    hold on;
                    scatter(spkpos(:,1),spkpos(:,2),'o','filled','CData',spk_pass_index,'SizeData',20);
                    line(...
                        [-floor(0.2*range(pos(:,1))/5)*5 0]+0.95*max(pos(:,1)),...
                        [1 1]*(min(pos(:,2)-0.1*range(pos(:,2)))),...
                        'Color','k','LineWidth',3);
                    text(0.95*max(pos(:,1)),...
                        min(pos(:,2))-0.15*range(pos(:,2)),...
                        [num2str(floor(0.2*range(pos(:,1))/5)*5) ' ' units],...
                        'VerticalAlignment','Cap','HorizontalAlignment','right');
                    axis off equal;
                    hold off;
                    
                case 3 % 3D
                    plot3(pos(:,1),pos(:,2),pos(:,3),'k','LineWidth',1)
                    spkpos = spk_pos(pos_ts,pos,spk_ts);
                    hold on;
                    scatter3(spkpos(:,1),spkpos(:,2),spkpos(:,3),...
                        'o','filled','CData',spk_pass_index);
                    l = floor(max(range(pos))*0.2/10)*10;
                    line([0 0 0;l 0 0]+min(pos(:,1))-0.01*range(pos(:,1)),...
                        [0 0 0;0 l 0]+min(pos(:,2))-0.01*range(pos(:,2)),...
                        [0 0 0;0 0 l]+min(pos(:,3))-0.01*range(pos(:,3)),...
                        'Color','k','LineWidth',3);
                    
                    text(min(pos(:,1))-0.02*range(pos(:,1)),...
                        min(pos(:,2))-0.02*range(pos(:,2)),...
                        min(pos(:,3))-0.02*range(pos(:,3)),...
                        [num2str(l) ' ' units],...
                        'VerticalAlignment','Cap','HorizontalAlignment','right');
                    hold off;
                    axis off equal;
                    view(3);
                otherwise % Can't plot
                    warning('pass_index:bad_dim_subplot','Cannot plot over 3d');
            end
        case 'Rate map'
            map =results.rate_map;
            switch size(pos,2)
                case 1 % 1D
                    bar(centers{:},map);
                    xlim(minmax(centers{:}));
                    ylim([0 max(map(:))*1.1]);
                    
                    temp = occupancy==0;
                    hold on;
                    while any(temp)
                        j = find(temp,1);
                        k = find(temp(j:end)==0,1)-2;
                        rectangle('Position',[centers{1}(j) 0 range(centers{1}(j:(j+k))) max(map(:))*1.1],'FaceColor',[1 1 1]*0.7,'LineStyle','none');
                        temp(j:(j+k))=0;
                    end
                    hold off;
                    ylabel(['State (' units ')']);
                    xlabel(['Rate (Hz)']);
                case 2 % 2D
                    pk = quantile(map(occupancy>0),0.99);
                    [cbar, clims] = smart_colorbar([0 pk], jet(255));
                    map(occupancy==0)=clims(1);
                    imagesc(centers{:},map');
                    colormap(cbar);
                    set(gca,'CLim',clims);
                    
                    hold on;
                    line(...
                        [-floor(0.2*range(pos(:,1))/5)*5 0]+0.95*max(pos(:,1)),...
                        [1 1]*(min(pos(:,2)-0.1*range(pos(:,2)))),...
                        'Color','k','LineWidth',3);
                    text(0.95*max(pos(:,1)),...
                        min(pos(:,2))-0.15*range(pos(:,2)),...
                        [num2str(floor(0.2*range(pos(:,1))/5)*5) ' ' units],...
                        'VerticalAlignment','cap','HorizontalAlignment','right');
                    text(max(pos(:,1)),...
                        max(pos(:,2)),...
                        [sprintf('%3.1f',pk) ' Hz'],...
                        'VerticalAlignment','bottom','HorizontalAlignment','right');
                    axis off equal;
                    set(gca,'YDir','normal');
                    hold off;
                    CMBHOME.Utils.freezeColors;
                case 3 % 3D
                    set(gcf,'Renderer','OpenGL');
                    pk = quantile(map(occupancy>0),0.99);
                    [cbar, clims] = smart_colorbar([0 pk], jet(255));
                    map(occupancy==0)=clims(1);
                    j = {0 1 2};
                    k = find(occupancy>0)';
                    view(3);
                    for k=k;
                        [j{1} j{2} j{3}] = ind2sub(size(map),k);
                        clr = round((map(k)-min(clims))/range(clims)*size(cbar,1))+1;
                        if clr>size(cbar,1), clr=size(cbar,1); end;
                        plotcube([1 1 1]*binside,...
                            cellfun(@(x,y)x(y),centers,j)-[0.5 0.5 0.5]*binside,...
                            cbar(clr,:),...
                            min((map(k)-min(clims))/range(clims),1)*0.5+0.02,...
                            'EdgeAlpha',0,...
                            'BackfaceCull',1);
                    end
                    l = floor(max(range(pos))*0.2/10)*10;
                    line([0 0 0;l 0 0]+min(pos(:,1))-0.01*range(pos(:,1)),...
                        [0 0 0;0 l 0]+min(pos(:,2))-0.01*range(pos(:,2)),...
                        [0 0 0;0 0 l]+min(pos(:,3))-0.01*range(pos(:,3)),...
                        'Color','k','LineWidth',3);
                    
                    text(min(pos(:,1))-0.02*range(pos(:,1)),...
                        min(pos(:,2))-0.02*range(pos(:,2)),...
                        min(pos(:,3))-0.02*range(pos(:,3)),...
                        [num2str(l) ' ' units],...
                        'VerticalAlignment','Cap','HorizontalAlignment','right');
                    
                    text(min(pos(:,1))-0.02*range(pos(:,1)),...
                        min(pos(:,2))-0.02*range(pos(:,2)),...
                        min(pos(:,3))-0.12*range(pos(:,3)),...
                        [sprintf('%3.1f',pk) ' Hz'],...
                        'VerticalAlignment','bottom','HorizontalAlignment','right');
                    hold off;
                    axis off equal;
                    map(occupancy==0)=clims(1);
                otherwise % Can't plot
                    warning('pass_index:bad_dim_subplot','Cannot plot over 3d');
            end
        case 'Field index map'
            if ~all(isnan(field_index_map(:)))
                switch size(pos,2)
                    case 1 %1D
                        bar(centers{:},field_index_map);
                        xlim(minmax(centers{:}));
                        ylim([0 max(field_index_map(:))*1.1]);
                        
                        temp = occupancy==0;
                        hold on;
                        while any(temp)
                            j = find(temp,1);
                            k = find(temp(j:end)==0,1)-2;
                            rectangle('Position',[centers{1}(j) 0 range(centers{1}(j:(j+k))) max(map(:))*1.1],'FaceColor',[1 1 1]*0.7,'LineStyle','none');
                            temp(j:(j+k))=0;
                        end
                        hold off;
                        ylabel(['State (' units ')']);
                        xlabel(['Field Index']);
                    case 2 %2D
                        [cbar, clims] = smart_colorbar([0 1], hot(255));
                        field_index_map(occupancy==0)=clims(1);
                        imagesc(centers{:},field_index_map');
                        colormap(cbar);
                        set(gca,'CLim',clims);
                        
                        hold on;
                        line(...
                            [-floor(0.2*range(pos(:,1))/5)*5 0]+0.95*max(pos(:,1)),...
                            [1 1]*(min(pos(:,2)-0.1*range(pos(:,2)))),...
                            'Color','k','LineWidth',3);
                        text(0.95*max(pos(:,1)),...
                            min(pos(:,2))-0.15*range(pos(:,2)),...
                            [num2str(floor(0.2*range(pos(:,1))/5)*5) ' ' units],...
                            'VerticalAlignment','cap','HorizontalAlignment','right');
                        axis off equal;
                        set(gca,'YDir','normal');
                        hold off;
                        CMBHOME.Utils.freezeColors;
                    case 3 %3D
                        set(gcf,'Renderer','OpenGL');
                        [cbar, clims] = smart_colorbar([0 1], hot(255));
                        field_index_map(occupancy==0)=clims(1);
                        j = {0 1 2};
                        k = find(occupancy>0)';
                        view(3);
                        for k=k;
                            [j{1} j{2} j{3}] = ind2sub(size(field_index_map),k);
                            clr = round((field_index_map(k)-min(clims))/range(clims)*size(cbar,1))+1;
                            if clr>size(cbar,1), clr=size(cbar,1); end;
                            plotcube([1 1 1]*binside,...
                                cellfun(@(x,y)x(y),centers,j)-[0.5 0.5 0.5]*binside,...
                                cbar(clr,:),...
                                min((field_index_map(k)-min(clims))/range(clims),1)*0.95+0.025,...
                                'EdgeAlpha',0,...
                                'BackfaceCull',1);
                        end
                        l = floor(max(range(pos))*0.2/10)*10;
                        line([0 0 0;l 0 0]+min(pos(:,1))-0.01*range(pos(:,1)),...
                            [0 0 0;0 l 0]+min(pos(:,2))-0.01*range(pos(:,2)),...
                            [0 0 0;0 0 l]+min(pos(:,3))-0.01*range(pos(:,3)),...
                            'Color','k','LineWidth',3);
                        
                        text(min(pos(:,1))-0.02*range(pos(:,1)),...
                            min(pos(:,2))-0.02*range(pos(:,2)),...
                            min(pos(:,3))-0.02*range(pos(:,3)),...
                            [num2str(l) ' ' units],...
                            'VerticalAlignment','Cap','HorizontalAlignment','right');
                        hold off;
                        axis off equal;
                    otherwise % Can't plot
                        warning('pass_index:bad_dim_subplot','Cannot plot over 3d');
                end
            end
        case 'Scatter plot'
            scatter(repmat(spk_pass_index,[2 1]),rad2deg([mod(spk_theta_phase,2*pi);mod(spk_theta_phase,2*pi)+2*pi]),'o','filled','SizeData',5);
            xlim([-1 1]);ylim([0 2*360]);
            xlabel('Pass Index');ylabel('LFP Phase (^o)');
            
            x = linspace(-1,1,750);
            phi = mod(2*pi*s*x+b,2*pi);
            k = find(abs(diff(phi))>pi);
            phi(k) = NaN;
            
            hold on;
            text(x(floor(750*0.75)),rad2deg(phi(floor(750*0.75))+2*pi)+25,['rho=' sprintf('%2.2f',rho)],'Color',[1 0 0],'FontWeight','bold','BackgroundColor','w');
            text(x(floor(750*0.75)),rad2deg(phi(floor(750*0.75))+2*pi)-25,['p=' sprintf('%2.2f',p)],'Color',[1 0 0],'FontWeight','bold','BackgroundColor','w');
            plot(x,rad2deg(phi+2*pi),'r','LineWidth',2);
            plot(x,rad2deg(phi),'r','LineWidth',2);
            
            hold off;
        case 'Density map'
            imagesc(dens_centers{1},rad2deg([dens_centers{2} dens_centers{2}+2*pi]),[density';density']);
            set(gca,'YDir','normal');
            xlabel('Pass Index');ylabel('LFP Phase (^o)');
            text(1,730,[sprintf('%3.2f',max(density(:))) 'Hz'],...
                'HorizontalAlignment','Right',...
                'VerticalAlignment','baseline');
            colormap jet;
            CMBHOME.Utils.freezeColors;
        otherwise
            warning('pass_index:bad_plot','Cannot recognize plot. Skipping...');
    end
    title(plots{i});
end


end

