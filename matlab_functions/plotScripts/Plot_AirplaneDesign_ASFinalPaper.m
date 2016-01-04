function Plot_AirplaneDesign_ASFinalPaper(PerfResults, DesignResults, environment, plane, params, flightdata, vars)
    %--------------------------------------------------------------------------
    %--- OVERVIEW PLOTS
    %--------------------------------------------------------------------------
    if(numel(vars(2).values)<2 ||numel(vars(1).values)<2) 
        display('ERROR: Not enough design variables specified, cannot plot results! Please specific at least the first two design variables.');
    else
        for i = 1:numel(vars(3).values)
            %Configure figure
            str = strcat(vars(3).shortname, '=', num2str(vars(3).values(i)));
            figure('Name',str);
            colormap(flipud(colormap(hot)));

            FontSize = 14;
            vmargin = 0.10;
            hmargin = 0.04;
            
            t_req = 6.9;
            
            %Convert data first
            for k = 1:numel(vars(2).values)
                for j = 1:numel(vars(1).values)
                    temp_t_exc(k,j) = PerfResults(i,k,j).t_excess;
                    temp_chargemargin(k,j) = PerfResults(i,k,j).t_chargemargin; 
                    temp_endurance(k,j) = PerfResults(i,k,j).t_endurance;
                    temp_m_total(k,j) = DesignResults(i,k,j).m_no_bat + DesignResults(i,k,j).m_bat;
                    if(isnan(temp_chargemargin(k,j))) temp_chargemargin(k,j) = 0.0; end %Remove NaNs to have smooth contour plotting below
                end
            end

            %NaN-check to avoid plotting errors below
            if(sum(~isnan(temp_endurance))==0) temp_endurance(:,:) = 0.0; end

            %PLOT : Excess time. Draws contour surface, and contour lines
            ax(1)=subplot_tight(2,2,1,[vmargin hmargin]);
            hold on
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_t_exc,120,'Linestyle','none');
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1]);
            clabel(c2,hc2,'LabelSpacing',144,'FontSize',FontSize);
            ylabel(vars(2).name);
            %xlabel(vars(1).name);
            title('Excess time T_{exc} [h]');
            cbar_handle(1) = colorbar;
            %set(get(cbar_handle,'ylabel'),'string','Excess time t_{exc} [h]')
            caxis([0,max(max(max( temp_t_exc(:,:))))])
            %set(gca,'DataAspectRatio',[1 4.5 1]);

            % PLOT : Charge margin
            ax(end+1)=subplot_tight(2,2,3,[vmargin hmargin]);
            hold on
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_chargemargin, 120,'Linestyle','none');
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1],'LineStyle','--');
            ylabel(vars(2).name);
            xlabel(vars(1).name);
            title('Charge margin T_{cm} [h]');
            cbar_handle(end+1)=colorbar;
            %set(get(cbar_handle,'ylabel'),'string','Charge margin t_{cm} [h]')
            caxis([0,max(max(max( temp_chargemargin(:,:))))])
            %set(gca,'DataAspectRatio',[1 4.5 1]);

            % PLOT : Endurance
            ax(end+1)=subplot_tight(2,2,2,[vmargin hmargin]);
            hold on
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_endurance, 120,'Linestyle','none');
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1],'LineStyle','--');
            %ylabel(vars(2).name);
            %xlabel(vars(1).name);
            title('Endurance T_{endur} [h]');
            cbar_handle(end+1)=colorbar;
            %set(get(cbar_handle,'ylabel'),'string','Endurance t_{endur} [h]')
            caxis([0,max(max(max( temp_endurance(:,:))))])
            %set(gca,'DataAspectRatio',[1 4.5 1]);

            % PLOT : Total airplane mass
            ax(end+1)=subplot_tight(2,2,4,[vmargin hmargin]);
            hold on
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_m_total, 120,'Linestyle','none');
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1],'LineStyle','--');
            %ylabel(vars(2).name);
            xlabel(vars(1).name);
            title('Total airplane mass m_{tot} [kg]');
            cbar_handle(end+1)=colorbar;
            %set(get(cbar_handle,'ylabel'),'string','Total mass m_{tot} [kg]')
            caxis([0,max(max(max( temp_m_total(:,:))))])
            %set(gca,'DataAspectRatio',[1 4.5 1]);
            
            %Set properties for all plots
            set(ax,'FontSize',FontSize);
            set(cbar_handle,'FontSize',FontSize);
            
            if(1)
                % Plot a cross to mark the chosen configuration?
                b_chosen = 5.6;
                m_bat_chosen = 2.9;
                for i = 1:numel(ax)
                    plot(ax(i),b_chosen,m_bat_chosen,'x','MarkerSize',9);
                end
            end
        end

        set(gcf, 'Color', 'w');
        %export_fig 3_excesstime.eps
    end    
    %--------------------------------------------------------------------------
    %--- SIMPLE TIME PLOTS
    %--------------------------------------------------------------------------
    for i = 1:numel(vars(3).values)
        for k = 1:numel(vars(2).values)
            for j = 1:numel(vars(1).values)
                if(abs(flightdata(i,k,j).b-5.6) <= eps(flightdata(i,k,j).b) && ...
                   abs(flightdata(i,k,j).m_bat-2.9) <= eps(flightdata(i,k,j).m_bat))
                    Plot_BasicSimulationTimePlot_ASJFR81hFlightPaper(flightdata(i,k,j), environment, params, plane);
                end
            end
        end
    end
    
    %export_fig test.eps;
end