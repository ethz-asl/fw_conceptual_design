function Plot_AirplaneDesign_ASFinalPaper_PlotOrderChanged(PerfResults, DesignResults, environment, plane, params, flightdata, vars)
    
    plotnr = 3;
    %set(0,'defaulttextinterpreter','latex');

    FontSize = 14;
    if(plotnr==3) Fontsize = 18;
    vmargin = 0.10;
    hmargin = 0.05;
            
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
            colormap(flipud(colormap(hot))); %colormap(jet);
            set(gcf, 'Color', 'w');
            
            %Convert data first
            for k = 1:numel(vars(2).values)
                for j = 1:numel(vars(1).values)
                    %display([i,' ',k,' ',j]);
                    temp_t_exc(k,j) = PerfResults(i,k,j).t_excess;
                    temp_min_SoC(k,j) = PerfResults(i,k,j).min_SoC;
                    temp_chargemargin(k,j) = PerfResults(i,k,j).t_chargemargin;
                    temp_endurance(k,j) = PerfResults(i,k,j).t_endurance;
                    temp_m_total(k,j) = 0.0;
                    if(~isempty(DesignResults)) DesignResults(i,k,j).m_no_bat + DesignResults(i,k,j).m_bat; end;
                    if(isnan(temp_chargemargin(k,j))) temp_chargemargin(k,j) = 0.0; end %Remove NaNs to have smooth contour plotting below
                end
            end

            %NaN-check to avoid plotting errors below
            if(sum(~isnan(temp_endurance))==0) temp_endurance(:,:) = 0.0; end
            %if(sum(~isnan(temp_chargemargin))==0) temp_chargemargin(:,:) = 0.0; end

            %PLOT1 : Excess time. Draws contour surface, and contour lines
            ax(1) = subplot_tight(2,2,1,[vmargin hmargin]);
            hold on
            ylabel_delta=0.0;
            if(plotnr==3) ylabel_delta=1.0; end;
            [c2,hc2] = contourf(vars(1).values,vars(2).values+ylabel_delta,temp_t_exc,[0.01:0.05:min(max(max(max(temp_t_exc))),24.0)],'Linestyle','none');
            [c2,hc2] = contour(vars(1).values,vars(2).values+ylabel_delta,temp_t_exc,'LevelStep',1.0,'LineColor',[0 0 0]);
            clabel(c2,hc2,'LabelSpacing',305,'FontSize',FontSize);
            title('Excess Time [h]');
            h = ylabel(vars(2).name);
            h = xlabel(vars(1).name);
            cbar_handle(1) = colorbar;
            caxis([0,min(max(max(max(temp_t_exc))),24.0)]);
            if(plotnr==3)
                xlabel('P_{in}/P_{in}^{nom} [-]');
                ylabel('P_{out}/P_{out}^{nom} [-]');
                caxis([0,10.10]);
                set(gca,'DataAspectRatio',[0.8 1 1]);
                set(get(cbar_handle(1),'ylabel'),'string','T_{exc} [h]')
            end
                
            % PLOT2 : Endurance
            ax(end+1) = subplot_tight(2,2,3,[vmargin hmargin]);
            hold on
            [c2,hc2] = contourf(vars(1).values,vars(2).values,temp_endurance, [0.1:0.1:max(max(max(temp_endurance)))],'Linestyle','none');
            [c2,hc2] = contour(vars(1).values,vars(2).values,temp_endurance,'LevelStep',2.0,'LineColor',[0 0 0],'ShowText','on');
            clabel(c2,hc2,'LabelSpacing',305,'FontSize',FontSize);
            title('Endurance [h]');
            ylabel(vars(2).name);
            xlabel(vars(1).name);
            cbar_handle(end+1) = colorbar;
            caxis([0,max(max(max(temp_endurance(:,:))))]);

            % PLOT3 : Charge margin
            ax(end+1) = subplot_tight(2,2,2,[vmargin hmargin]);
            hold on
            [c2,hc2] = contourf(vars(1).values,vars(2).values,temp_chargemargin, [0.01:0.05:min(max(max(max(temp_chargemargin))),24.0)],'Linestyle','none');
            [c2,hc2] = contour(vars(1).values,vars(2).values,temp_chargemargin,'LevelStep',1.0,'LineColor',[0 0 0],'ShowText','on');
            title('Charge margin [h]');
            ylabel(vars(2).name);
            xlabel(vars(1).name);
            cbar_handle(end+1) = colorbar;
            caxis([0,min(max(max(max(temp_chargemargin))),24.0)]);            

            %PLOT4: Optional plots can be chosen here

        %     PLOT4.1: Structural Mass
        %     subplot(2,2,4)
        %     contourf(designVars.b_array,...
        %         designVars.m_bat_array, m_struct(:,:,l)',200,...
        %         'Linestyle','none');
        %     xlabel('Wingspan (b) [m]')
        %     ylabel('Battery mass (m\_bat) [kg]');
        %     title(['Structural Mass [kg] (AR=' num2str(AR) ')']);
        %     caxis([min(min( m_struct(:,:,l))),max(max(m_struct(:,:,l)))])
        %     colorbar

    %         % PLOT4.1: Combined Performace (T_exc+T_ChargeMargin)
    %         subplot(2,2,4)
    %         a=0.6;
    %         b=1.0-a;
    %         comb_perf=a*results.t_chargemargin+b*results.t_excess;
    %         for i=1:size(results.flightdata_array,1)
    %             for j=1:size(results.flightdata_array,2)
    %                 for k=1:size(results.flightdata_array,3)
    %                     if results.flightdata_array(i,j,k).minSoC<0.1
    %                         comb_perf(i,j,k)=NaN;
    %                     end
    %                 end
    %             end
    %         end
    % 
    %         contourf(designVars.b_array,...
    %             designVars.m_bat_array, comb_perf(:,:,l)',200,...
    %             'Linestyle','none');
    %         xlabel('Wingspan (b) [m]')
    %         ylabel('Battery mass (m\_bat) [kg]');
    %         title(['a*{t_cm}+b*t_{exc} [kg] (AR=' num2str(AR) ')']);
    %         caxis([min(min( comb_perf(:,:,l))),max(max(comb_perf(:,:,l)))])
    %         colorbar

            % PLOT4.3: Total airplane mass
            ax(end+1) = subplot_tight(2,2,4,[vmargin hmargin]);
            hold on
            [c2,hc2] = contourf(vars(1).values,vars(2).values,temp_m_total, [0.01:0.01:max(max(max(temp_m_total)))],'Linestyle','none');
            [c2,hc2] = contour(vars(1).values,vars(2).values,temp_m_total,'LevelStep',1.0,'LineColor',[0 0 0],'ShowText','on');
            title('Total mass [kg]');
            ylabel(vars(2).name);
            xlabel(vars(1).name);
            cbar_handle(end+1) = colorbar;
            caxis([0,max(max(max(temp_m_total(:,:))))]);
            
            %Set properties for all plots
            set(ax,'FontSize',FontSize);
            set(cbar_handle,'FontSize',FontSize);
        end
    end

    %--------------------------------------------------------------------------
    %--- DETAILED TIME PLOTS
    %--------------------------------------------------------------------------
    %DetailedAnalysisForm({designVars.AR_array,designVars.b_array,designVars.m_bat_array,results.flightdata_array,results.t_excess,results.t_endurance,results.t_chargemargin});

    %--------------------------------------------------------------------------
    %--- SIMPLE TIME PLOTS
    %--------------------------------------------------------------------------
    % Allows to choose one configuration and plot detailed data.
    for i = 1:numel(vars(3).values)
        for k = 1:numel(vars(2).values)
            for j = 1:numel(vars(1).values)
                if(flightdata(i,k,j).b == 4.5 && flightdata(i,k,j).m_bat == 7.0)
                    Plot_BasicSimulationTimePlot(flightdata(i,k,j), environment, params, plane);
                end
            end
        end
    end
    
end     