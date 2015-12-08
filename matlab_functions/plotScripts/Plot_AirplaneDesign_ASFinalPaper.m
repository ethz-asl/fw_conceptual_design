function Plot_AirplaneDesign_ASFinalPaper(PerfResults, DesignResults, environment, plane, parameters, flightdata, vars)
    %--------------------------------------------------------------------------
    %--- OVERVIEW PLOTS
    %--------------------------------------------------------------------------
    if(numel(vars(2).values)<2 ||numel(vars(1).values)<2) 
        display('ERROR: Not enough design variables specified, cannot plot results! Please specific at least the first two design variables.');
    else
        for i = 1:numel(vars(3).values)
            %Set figure title
            str = strcat(vars(3).shortname, '=', num2str(vars(3).values(i)));
            figure('Name',str);

            %Convert data first
            for k = 1:numel(vars(2).values)
                for j = 1:numel(vars(1).values)
                    display([i,' ',k,' ',j]);
                    temp_t_exc(k,j) = PerfResults(i,k,j).t_excess;
                    temp_chargemargin(k,j) = PerfResults(i,k,j).t_chargemargin; 
                    temp_endurance(k,j) = PerfResults(i,k,j).t_endurance;
                    temp_m_total(k,j) = DesignResults(i,k,j).m_no_bat + DesignResults(i,k,j).m_bat;
                    if(isnan(temp_chargemargin(k,j))) temp_chargemargin(k,j) = 0.0; end %Remove NaNs to have smooth contour plotting below
                    if(isnan(temp_endurance(k,j))) temp_endurance(k,j) = 0.0; end %Remove NaNs to have smooth contour plotting below
                end
            end

            %NaN-check to avoid plotting errors below
           % if(isnan(temp_chargemargin)==0) temp_chargemargin(:,:) = 0.0; end

            t_req = 6.9;
            FontSize = 14;
            %PLOT : Excess time. Draws contour surface, and contour lines
            subplot_tight(2,2,1,[0.06 0.06])
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_t_exc,120,'Linestyle','none');
            hold on
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1]);
            clabel(c2,hc2,'LabelSpacing',144,'FontSize',FontSize);
            ylabel(vars(2).name);
            %xlabel(vars(1).name);
            cbar_handle = colorbar;
            set(get(cbar_handle,'ylabel'),'string','Excess time t_{exc} [h]')
            caxis([0,max(max(max( temp_t_exc(:,:))))])
            set(gca,'DataAspectRatio',[1 4.5 1]);
            set(gca,'FontSize',FontSize);
            set(cbar_handle,'FontSize',FontSize);
            colormap(flipud(colormap(hot)));

            % PLOT : Charge margin
            subplot_tight(2,2,3,[0.06 0.06])
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_chargemargin, 120,'Linestyle','none');
            hold on
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1],'LineStyle','--');
            ylabel(vars(2).name);
            xlabel(vars(1).name);
            cbar_handle=colorbar;
            set(get(cbar_handle,'ylabel'),'string','Charge margin t_{cm} [h]')
            caxis([0,max(max(max( temp_chargemargin(:,:))))])
            set(gca,'DataAspectRatio',[1 4.5 1]);
            set(gca,'FontSize',FontSize);
            set(cbar_handle,'FontSize',FontSize);
            colormap(flipud(colormap(hot)));

            % PLOT : Endurance
            subplot_tight(2,2,2,[0.06 0.06])
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_endurance, 120,'Linestyle','none');
            hold on
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1],'LineStyle','--');
            ylabel(vars(2).name);
            %xlabel(vars(1).name);
            cbar_handle=colorbar;
            set(get(cbar_handle,'ylabel'),'string','Endurance t_{endur} [h]')
            caxis([0,max(max(max( temp_endurance(:,:))))])
            set(gca,'DataAspectRatio',[1 4.5 1]);
            set(gca,'FontSize',FontSize);
            set(cbar_handle,'FontSize',FontSize);
            colormap(flipud(colormap(hot)));
            ylabel(vars(2).name);
            xlabel(vars(1).name);

            % PLOT : Total airplane mass
            subplot_tight(2,2,4,[0.06 0.06])
            [c2,hc2]=contourf(vars(1).values,vars(2).values,temp_m_total, 120,'Linestyle','none');
            hold on
            [c2,hc2]=contour(vars(1).values,vars(2).values,temp_t_exc,[t_req t_req],'LineColor',[0 0 1],'LineStyle','--');
            ylabel(vars(2).name);
            xlabel(vars(1).name);
            cbar_handle=colorbar;
            set(get(cbar_handle,'ylabel'),'string','Total mass m_{tot} [kg]')
            caxis([0,max(max(max( temp_m_total(:,:))))])
            set(gca,'DataAspectRatio',[1 4.5 1]);
            set(gca,'FontSize',FontSize);
            set(cbar_handle,'FontSize',FontSize);
            colormap(flipud(colormap(hot)));
            ylabel(vars(2).name);
            xlabel(vars(1).name);
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
                if(flightdata(i,k,j).b == 4.5 && flightdata(i,k,j).m_bat == 7.0)
                    Plot_BasicSimulationTimePlot(flightdata(i,k,j), environment, parameters, plane);
                end
            end
        end
    end
    
    %export_fig test.eps;
end