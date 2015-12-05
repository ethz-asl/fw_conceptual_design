function Plot_DesignSpaceEnvironment_ASFinalPaper
    %Plotting

    n=0;

    figure
    texc_reduced(:,:)=t_excess(:,1,:);
    [c2,hc2]=contourf(clearness_array,turbulence_array+1,texc_reduced(:,:),...
        200,'Linestyle','none','ShowText','off');
    colormap(hot); %Set an inverted version of the "hot" colormap
    cmap = colormap;
    cmap = flipud(cmap);
    colormap(cmap);
    xlabel('P_{in}/P_{in,nom} [-]');
    ylabel('P_{out}/P_{out,nom} [-]');
    caxis([0,9]);
    set(gca,'DataAspectRatio',[0.8 1 1]);
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [0 0 300 300]);
    cbar_handle=colorbar;
    xbounds = xlim;
    set(gca,'XTick',fliplr(xbounds(2):-0.1:xbounds(1)));
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):0.1:ybounds(2));
    set(get(cbar_handle,'ylabel'),'string','t_{exc} [h]')
    if(environment.month==6)
        export_fig 5a_texc.eps
    else
        export_fig 5c_texc.eps
    end

    figure
    texc_reduced(:,:)=t_excess(:,2,:);
    [c2,hc2]=contourf(clearness_array,turbulence_array+1,texc_reduced(:,:),...
        200,'Linestyle','none','ShowText','off');
    colormap(hot); %Set an inverted version of the "hot" colormap
    cmap = colormap;
    cmap = flipud(cmap);
    colormap(cmap);
    xlabel('P_{in}/P_{in,nom} [-]');
    ylabel('P_{out}/P_{out,nom} [-]');
    caxis([0,9]);
    set(gca,'DataAspectRatio',[0.8 1 1]);
    set(gcf, 'Color', 'w');
    set(gcf, 'Position', [100 100 300 300]);
    cbar_handle=colorbar;
    set(get(cbar_handle,'ylabel'),'string','t_{exc} [h]')
    xbounds = xlim;
    ybounds = ylim;
    set(gca,'YTick',ybounds(1):0.1:ybounds(2));
    set(gca,'XTick',fliplr(xbounds(2):-0.1:xbounds(1)));
    if(environment.month==6)
        export_fig 5b_texc.eps
    else
        export_fig 5d_texc.eps
    end

    % for m_bat=m_bat_min:m_bat_step:m_bat_max
    %     if(length(clearness_array)<2 ||length(turbulence_array)<2) continue; end
    %     
    %     n=n+1;
    %     
    %     subplot_tight(2,1,n,[0.08 0.1])
    %     texc_reduced(:,:)=t_excess(:,n,:);
    %     [c2,hc2]=contourf(clearness_array,turbulence_array+1,texc_reduced(:,:),...
    %         200,'Linestyle','none','ShowText','off');
    %     colormap(hot); %Set an inverted version of the "hot" colormap
    %     cmap = colormap;
    %     cmap = flipud(cmap);
    %     colormap(cmap);
    %     ylabel('P_{out}/P_{out,nom} [-]')
    %     xlabel('P_{in}/P_{in,nom} [-]');
    %     %title(['Excess Time [h] (m_{bat}=' num2str(m_bat) 'kg)']);
    %     %caxis([-15,15])
    %     caxis([0,max(max(max( t_excess(:,:,:))))])
    %     cbar_handle=colorbar
    %     set(get(cbar_handle,'ylabel'),'string','t_{exc} [h]')
    %     
    % %     nc = get(hc2,'Children');
    % %     temp = 100;
    % %     select_i=zeros(length(nc),1);
    % %     for i = 1:length(nc)
    % %        ud1 = get(nc(i),'UserData');   
    % %        if (abs(ud1) < temp)
    % %            temp = abs(ud1);
    % %        end
    % %     end
    % %     filler=1;
    % %     for i = 1:length(nc)
    % %        ud1 = get(nc(i),'UserData');   
    % %        if (abs(ud1) == temp)
    % %            select_i(filler)=i;
    % %            filler=filler+1;
    % %        end
    % %     end
    % %     j=1;
    % %     if i>0
    % %         while(select_i(j)~=0)
    % %             set(nc(select_i(j)),'Linestyle','-');
    % %             set(nc(select_i(j)),'LineWidth',2);
    % %             j=j+1;
    % %         end
    % %     end
    % 
    % 
    % end
end