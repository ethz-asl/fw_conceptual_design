function Plot_DesignSpace_ASFinalPaper

    %--------------------------------------------------------------------------
    %--- OVERVIEW PLOTS
    %--------------------------------------------------------------------------
    l=1;
    i=0;

    for AR=AR_min:AR_step:AR_max
        if(length(b_array)<2 ||length(m_bat_array)<2) continue; end

        i=i+1;
        h12(l)=figure(l+1); %TOCHANGE

        %Draw the contour surface
        %subplot(2,2,1)
        [c2,hc2]=contourf(b_min:b_step:b_max,...
            m_bat_min:m_bat_step:m_bat_max, t_excess(:,:,l)',...
            200,'Linestyle','none');
        colormap(flipud(colormap(hot)));
        %xlabel('Wingspan b [m]')
        ylabel('Battery mass m_{bat} [kg]');
        %title(['Excess time t_{exc} [h] (\lambda=' num2str(AR) ')']);
        %caxis([-15,15])
        caxis([0,max(max(max( t_excess(:,:,:))))])
        set(gca,'DataAspectRatio',[1 4.5 1]);
        cbar_handle=colorbar
        set(get(cbar_handle,'ylabel'),'string','Excess time t_{exc} [h]')

        %Draw the contour lines at t_exc=0
    %     nc = get(hc2,'Children');
    %     temp = 100;
    %     select_i=zeros(length(nc),1);
    %     for i = 1:length(nc)
    %        ud1 = get(nc(i),'UserData');   
    %        if (abs(ud1) < temp)
    %            temp = abs(ud1);
    %        end
    %     end
    %     filler=1;
    %     for i = 1:length(nc)
    %        ud1 = get(nc(i),'UserData');   
    %        if (abs(ud1) == temp)
    %            select_i(filler)=i;
    %            filler=filler+1;
    %        end
    %     end
    %     j=1;
    %     if i>0
    %         while(select_i(j)~=0)
    %             set(nc(select_i(j)),'Linestyle','-');
    %             set(nc(select_i(j)),'LineWidth',1);
    %             set(nc(select_i(j)),'Color',white);
    %             j=j+1;
    %         end
    %     end

        %Draw the contour lines at t_exc=t_req
        nc = get(hc2,'Children');
        t_req = 6.9;
        select_i=zeros(length(nc),1);
        for i = 1:length(nc)
           ud1 = get(nc(i),'UserData');   
           if (ud1 < t_req)
               t_req = ud1;
               break;
           end
        end
        filler=1;
        for i = 1:length(nc)
           ud1 = get(nc(i),'UserData');   
           if (abs(ud1) == t_req)
               select_i(filler)=i;
               filler=filler+1;
           end
        end
        j=1;
        if i>0
            while(select_i(j)~=0)
                %set(nc(select_i(j)),'Linestyle','-');
                %set(nc(select_i(j)),'LineWidth',2);
                %set(nc(select_i(j)),'Color','blue');
                coords(j).x=get(nc(select_i(j)),'XData'); %Get coordinates of the contour
                coords(j).y=get(nc(select_i(j)),'YData'); %Get coordinates of the contour
                line(coords(j).x,coords(j).y,'Color','blue','LineStyle','-','LineWidth',1);
                j=j+1;
            end
        end

        set(gcf, 'Color', 'w');
        export_fig 3_excesstime.eps

    %     figure
    %     subplot(2,2,2)
    %     contourf(b_min:b_step:b_max,...
    %         m_bat_min:m_bat_step:m_bat_max, t_endurance(:,:,l)',...
    %         200,'Linestyle','none');
    %     xlabel('Wingspan (b) [m]')
    %     ylabel('Battery mass (m\_bat) [kg]');
    %     title(['Endurance [h] (AR=' num2str(AR) ')']);
    %     caxis([0,48])
    %     colorbar

        %subplot(2,2,3)
        figure
        contourf(b_min:b_step:b_max,...
            m_bat_min:m_bat_step:m_bat_max, t_chargemargin(:,:,l)',...
            200,'Linestyle','none');
        colormap(flipud(colormap(hot)));

        %Copy the  contour from the t_excess plot here
        for j=1:size(coords,1)
            line(coords(j).x,coords(j).y,'Color','blue','LineStyle','--');
        end

        xlabel('Wingspan b [m]')
        ylabel('Battery mass m_{bat} [kg]');    
        %title(['Charge margin t_{cm} [h] (\lambda=' num2str(AR) ')']);
        caxis([0,max(max(max( t_chargemargin(:,:,:))))])
        set(gca,'DataAspectRatio',[1 4.5 1]);
        cbar_handle=colorbar
        set(get(cbar_handle,'ylabel'),'string','Charge margin t_{cm} [h]')
        set(gcf, 'Color', 'w');
        export_fig 4_chargemargin.eps

    % %   Structural Mass
    %     subplot(2,2,4)
    %     contourf(b_min:b_step:b_max,...
    %         m_bat_min:m_bat_step:m_bat_max, m_struct(:,:,l)',200,...
    %         'Linestyle','none');
    %     xlabel('Wingspan (b) [m]')
    %     ylabel('Battery mass (m\_bat) [kg]');
    %     title(['Structural Mass [kg] (AR=' num2str(AR) ')']);
    %     caxis([min(min( m_struct(:,:,l))),max(max(m_struct(:,:,l)))])
    %     colorbar

    %   Combined Performance: T_exc+T_ChargeMargin
    %     subplot(2,2,4)
    %     a=0.6;
    %     b=1.0-a;
    %     comb_perf=a*t_chargemargin+b*t_excess;
    %     for i=1:size(flightdata_array,1)
    %         for j=1:size(flightdata_array,2)
    %             for k=1:size(flightdata_array,3)
    %                 if flightdata_array(i,j,k).minSoC<0.1
    %                     comb_perf(i,j,k)=NaN;
    %                 end
    %             end
    %         end
    %     end
    %     
    %     contourf(b_min:b_step:b_max,...
    %         m_bat_min:m_bat_step:m_bat_max, comb_perf(:,:,l)',200,...
    %         'Linestyle','none');
    %     xlabel('Wingspan (b) [m]')
    %     ylabel('Battery mass (m\_bat) [kg]');
    %     title(['a*{t_cm}+b*t_{exc} [kg] (AR=' num2str(AR) ')']);
    %     caxis([min(min( comb_perf(:,:,l))),max(max(comb_perf(:,:,l)))])
    %     colorbar

        l=l+1;
    end

    %--------------------------------------------------------------------------
    %--- DETAILED TIME PLOTS
    %--------------------------------------------------------------------------
    %DetailedAnalysisForm({AR_array,b_array,m_bat_array,flightdata_array,t_excess,t_endurance,t_chargemargin});

    %--------------------------------------------------------------------------
    %--- SIMPLE TIME PLOTS
    %--------------------------------------------------------------------------
    figure
    [ax]=plotyy([flightdata_array(1,1,1).t_array'/3600,flightdata_array(1,1,1).t_array'/3600],...
            [flightdata_array(1,1,1).P_solar_array',flightdata_array(1,1,1).P_elec_tot_array'],...
           flightdata_array(1,1,1).t_array/3600,...
            flightdata_array(1,1,1).bat_array/3600);
    l1=legend('P_{Solar}[W]','P_{elec,tot}[W]','E_{Bat} [Wh]');
    xlabel('Time [h]');
    ylabel(ax(1),'Power[W]');
    ylabel(ax(2),'Energy[Wh]');
    ylim(ax(2),[0 max(flightdata_array(1,1,1).bat_array/3600)*1.1]);
    %set(gcf, 'Position', [100 100 150 150]);
    saveas(gcf, 'test.eps');

    %end use double core- processing
     %matlabpool close
     %print(h12(1),'test2.emf','-dmeta','-painters')
 
end